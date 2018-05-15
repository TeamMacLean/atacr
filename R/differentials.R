#' gets t-statistic for two vectors of data, x and y
#' @param data matrix of sample data
#' @param indices indices selected by boot::boot
#' @return t the t statistic from Student's t-test or NA if error
get_t <- function(data,indices){
  d <- data[indices]

  e <- length(d)
  f <- floor(e/2)
  x <- d[1:f]
  y <- d[(f+1):e]


  stat <- tryCatch({
    t.test(x,y)$statistic
  },
    warning = function(w){
      return(NA)
    },
    error = function(e){
      return(NA)
    },
    finally = {}
  )
  return(stat)
}


#'runs bootstrap t test, wrapper required for boot::boot function
#' @param data matrix of sample data
#' @param iterations number of bootstrap iterations to run
#' @return vector of 2 columns, observed value t statisitc and p, calculated as proportion of bootstrap iterations greater than original t
bootstrap_t <- function(data, iterations=10){
  boot_res <- boot::boot(data, statistic = get_t, R = iterations)
  original <- boot_res$t0
  bootstraps <- boot_res$t
  p <- (sum(bootstraps > original) / iterations)
  if ( is.nan(original) | is.na(original) ) {
    p <- original
  }
  else if (original < 0) {
    p <- sum(bootstraps < original) / iterations
  }
  return(c(original, p))
}

select_comparisons <- function(data, treatment_a, treatment_b, which = "bait_windows"){
    l <- list()
    sample_matrix <- SummarizedExperiment::assay(data[[which]])
    treatment_a_cols <- data$sample_names[which(data$treatments == treatment_a) ]
    treatment_b_cols <- data$sample_names[which(data$treatments == treatment_b) ]
    l$treatment_a_data <- sample_matrix[,treatment_a_cols]
    l$treatment_b_data <- sample_matrix[,treatment_b_cols]
    return(l)
}

get_means <- function(data){

  mean_count_a <- apply(data$comparisons$treatment_a_data, 1, mean)
  mean_count_b <- apply(data$comparisons$treatment_b_data, 1, mean)
  result <- data.frame(
    c1 = mean_count_a,
    c2 = mean_count_b
  )
  colnames(result) <- c(paste0("mean_", data$treatment_a_name), paste0("mean_", data$treatment_b_name) )
  return( result )
}

get_sd <- function(data){
  sd_a <- apply(data$comparisons$treatment_a_data, 1, sd)
  sd_b <- apply(data$comparisons$treatment_b_data, 1, sd)
  result <- data.frame(
    c1 = sd_a,
    c2 = sd_b
  )
  colnames(result) <- c(paste0("sd_", data$treatment_a_name), paste0("sd_", data$treatment_b_name) )
  return( result )
}

get_fc <- function(data){
  means <- get_means(data)
  return(data.frame(log2_fold_change = log2(means[,1] / means[,2])))
}

#' selects appropriate columns and names from a
#' @param data an atacr object
#' @param treatment_a string naming the first treatment (numerator)
#' @param treatment_b string naming the second treatment (denominator)
#' @param which subset to work on Default = NULL
#' @return list of data to be calculated with
select_data <- function(data, treatment_a, treatment_b, which = NULL){

  comparison_list <- select_comparisons(data, treatment_a, treatment_b, which = which)
  comparison_matrix <- cbind(comparison_list$treatment_a, comparison_list$treatment_b )

  return(
    list(
      counts = comparison_matrix,
      comparisons = comparison_list,
      treatment_a_names = data$sample_names[which(data$treatments == treatment_a)],
      treatment_b_names = data$sample_names[which(data$treatments == treatment_b)],
      treatment_a_name = treatment_a,
      treatment_b_name = treatment_b
      )
    )

}

check_data <- function(d, treatment_a, treatment_b){

  if( length(d$treatment_a_names) < 3 | length(d$treatment_b_names) < 3  ){
    message <- paste("Need at least 3 replicates to perform estimate bootstrap t value. Have", length(d$treatment_a_names), "for", treatment_a, "and", length(d$treatment_b_names), "for", treatment_b)
    stop(message)
  }

  if(length(d$treatment_a_names) != length(d$treatment_b_names) ){
    message <- paste("Must have equal number of replicates in each treatment. Have", length(d$treatment_a_names), "for", treatment_a, "and", length(d$treatment_b_names), "for", treatment_b)
    stop(message)
  }

}

#' Estimate FDR and significantly different windows
#' @export
#' @param data an atacr object
#' @param treatment_a the first treatment to consider
#' @param treatment_b the second treatment to consider
#' @param which the subset of windows to consider
#' @param iterations the number of bootstrap iterations to perform
#' @param fdr_level the level at which to mark FDR as significant
#' @return dataframe of counts and statistics
estimate_fdr <- function(data, treatment_a, treatment_b, which = "bait_windows", iterations=10,fdr_level=0.05){


  d <- select_data(data, treatment_a, treatment_b, which)
  check_data(d, treatment_a, treatment_b)

  working_df <- as.data.frame(d$counts)
  row.names(working_df) <- rownames(d$counts)

  selected_df <- working_df[rowSums(working_df) > 0,]

  selected_result <- apply(selected_df, 1, bootstrap_t, iterations = iterations)
  #colnames(selected_result) <- c("t", "fdr")
  #selected_result <- apply(selected_df, 1, bayes_t, treatment_a_names = d$treatment_a_names, treatment_b_names = d$treatment_b_names)

 selected_result <- as.data.frame(t(selected_result)) %>%
  dplyr::rename("fdr" = V2) %>%
    dplyr::mutate(window = colnames(selected_result))



  working_df$window <- row.names(working_df)

  result <- dplyr::left_join(working_df, selected_result, by = "window") %>%
    dplyr::mutate(is_sig = fdr <= fdr_level) %>%
    dplyr::bind_cols( get_means(d) ) %>%
    dplyr::bind_cols( get_sd(d) ) %>%
    dplyr::bind_cols( get_fc(d))
  return(result)

}

#' Estimate FDR and significantly different windows for many experiments
#' @export
#' @param data an atacr object
#' @param common_control the treatment to consider the control for all other treatments
#' @param which the subset of windows to consider
#' @param iterations the number of bootstrap iterations to perform
#' @param fdr_level the level at which to mark FDR as significant
estimate_fdr_multiclass <- function(data, common_control, which = "bait_windows", iterations = 10,fdr_level = 0.05) {
  treatments <- data$treatments[data$treatments != common_control]
  control <- rep(common_control, length(treatments))
  comparisons <- cbind(treatments, control)

  r <- list()
  for (i in 1:nrow(comparisons)) {
    tr <- comparisons[i,][1]
    ct <- comparisons[i,][2]

    df <- atacr::estimate_fdr(data,
                    tr,
                    ct,
                    which = which,
                    iterations = iterations,
                    fdr_level = fdr_level)
         df$a <- rep(tr, nrow(df))
         df$b <- rep(ct, nrow(df))
         colnames(df)[grep("mean_", colnames(df))] <- c("mean_a", "mean_b")
         colnames(df)[grep("sd_", colnames(df))] <- c("sd_a", "sd_b")
         df <- df[, c("window",  "fdr", "is_sig", "mean_a", "mean_b", "sd_a", "sd_b", "log2_fold_change", "a", "b")]
    r[[i]] <- df
    }
  return(do.call(rbind, r))

}

bayes_t <- function(counts, treatment_a_names, treatment_b_names){

  a <- counts[treatment_a_names]
  b <- counts[treatment_b_names]
  bf <- BayesFactor::ttestBF(a,b)
  return(bf@bayesFactor$bf)
}


#' Estimate Bayes Factor and significantly different windows
#' @export
#' @param atacr an atacr object
#' @param treatment_a the first treatment to consider
#' @param treatment_b the second treatment to consider
#' @param which the subset of windows to consider
#' @param factor the BayesFactor at which to mark window as significant
#' @return a dataframe of counts and statistics
estimate_bayes_factor <- function(atacr, treatment_a, treatment_b, which = "bait_windows", factor = 4){

  d <- select_data(atacr, treatment_a, treatment_b, which)
  check_data(d, treatment_a, treatment_b)

  working_df <- as.data.frame(d$counts)
  row.names(working_df) <- rownames(d$counts)

  selected_df <- working_df[rowSums(working_df) > 0,]

  selected_result <- apply(selected_df, 1, bayes_t, treatment_a_names = d$treatment_a_names, treatment_b_names = d$treatment_b_names)

   selected_result <- data.frame(
     bayes_factor = selected_result,
     window = row.names(selected_df)
   )

  working_df$window <- row.names(working_df)
   result <- dplyr::left_join(working_df, selected_result, by = "window") %>%
     dplyr::mutate(is_sig = bayes_factor >= factor) %>%
     dplyr::bind_cols( get_means(d) ) %>%
     dplyr::bind_cols( get_sd(d) ) %>%
     dplyr::bind_cols( get_fc(d))
  return(result)
}

#' Estimate BayesFactor and mark significantly different windows for many experiments
#' @export
#' @param data an atacr object
#' @param common_control the treatment to consider the control for all other treatments
#' @param which the subset of windows to consider
#' @param factor the BayesFactor to consider significant
#' @return a dataframe of counts and statistics
estimate_bayes_factor_multiclass <- function(data, common_control, which = "bait_windows", factor = 4) {
  treatments <- unique(data$treatments[data$treatments != common_control])
  control <- rep(common_control, length(treatments))
  comparisons <- cbind(treatments, control)
  r <- list()
  for (i in 1:nrow(comparisons)) {
    tr <- comparisons[i,][1]
    ct <- comparisons[i,][2]
    df <- estimate_bayes_factor(data,
      tr,
      ct,
      which = which,
      factor = factor)
    df$a <- rep(tr, nrow(df))
    df$b <- rep(ct, nrow(df))
    colnames(df)[grep("mean_", colnames(df))] <- c("mean_a", "mean_b")
    colnames(df)[grep("sd_", colnames(df))] <- c("sd_a", "sd_b")
    df <- df[, c("window",  "bayes_factor", "is_sig", "mean_a", "mean_b", "sd_a", "sd_b", "log2_fold_change", "a", "b")]
    r[[i]] <- df
  }

  return(do.call(rbind, r))

}
#' Estimate differential window counts  and mark significantly different windows using edgeR exact method for two samples
#' @export
#' @param atacr an atacr object
#' @param common_control the treatment to consider the control for all other treatments
#' @param which the subset of windows to consider
#' @param sig_level the p_value to consider significant
#' @return a dataframe of counts and statistics
edgeR_exact <- function(atacr, which = "bait_windows", treatment_a = NULL, treatment_b = NULL, remove_zeros = FALSE, sig_level = 0.05 ){

  data <- select_data(atacr, treatment_a, treatment_b, which)
  working_df <- as.data.frame(data$counts)
  row.names(working_df) <- rownames(data$counts)

  group <- c(rep(treatment_a, length(data$treatment_a_names)), rep(treatment_b, length(data$treatment_b_names)) )

  dg <- edgeR::DGEList(data$counts, group = group, remove.zeros = remove_zeros)
  dg <- edgeR::estimateDisp(dg)
  et <- edgeR::exactTest(dg)
  names <- rownames(et$table)

  selected_result <- data.frame(
    window = rownames(et$table),
    p_value = et$table$PValue
  )

  working_df$window <- row.names(working_df)

  result <- dplyr::left_join(working_df, selected_result, by = "window") %>%
    dplyr::mutate(is_sig = p_value <= sig_level) %>%
    dplyr::bind_cols( get_means(data) ) %>%
    dplyr::bind_cols( get_sd(data) ) %>%
    dplyr::bind_cols( get_fc(data))
  return(result)
}
#' Estimate differential window counts  and mark significantly different windows using edgeR glmFIT method for multiple samples with common control
#' @export
#' @param data an atacr object
#' @param treatment_a the first treatment to consider
#' @param treatment_b the second treatment to consider
#' @param which the subset of windows to consider
#' @param remove_zeros apply edgeR remove.zeros argument
#' @return a list of "DGELRT" objects for each comparison
edgeR_multiclass <- function(data, common_control, which = "bait_windows", sig_level = 0.05, remove_zeros = FALSE){

  ctrl_idcs <- which(data$treatments == common_control)
  other_idcs <- which(data$treatments != common_control)
  new_order <- c(ctrl_idcs, other_idcs)

  treatments <- as.factor(data$treatments[ new_order ])
  samples <- data$sample_names[ new_order ]


  df <- data.frame(sample = samples, treatment = as.factor(as.numeric(treatments)))
  design <- model.matrix(~treatment, data = df)
  num_levels <- nlevels(as.factor(unique(treatments)))

  dglist <- edgeR::DGEList(SummarizedExperiment::assay(data[[which]]), remove.zeros = remove_zeros)

  dglist <- edgeR::estimateDisp(dglist, design)
  fit <- edgeR::glmQLFit(dglist, design)

  dgelrts <- list()

  for (i in 2:num_levels) {
    curr_t <- unique(data$treatments[ new_order ])[i]
    dgelrts[[curr_t]] <- edgeR::glmQLFTest(fit, coef = i)
  }

  return(dgelrts)
  # result <- list()
  #
  # for(n in names(dgelrts)){
  #   tb <- dgelrts[[n]]$table
  #   df <-  data.frame(
  #     window = rownames(tb),
  #     p_value = tb$PValue,
  #     f = tb$F
  #   )
  #
  #
  #   dlist <- select_data(data, n, common_control, which)
  #   df$is_sig <- (df$p_value <= sig_level)
  #
  #   df <- cbind(df, get_means(dlist$comparisons))
  #
  #   #add sd
  #   df <- cbind(df, get_sd(dlist$comparisons))
  #   #add log2 fc
  #   df <- get_fc(df)
  #   df$a <- rep(n, nrow(df))
  #   df$b <- rep(common_control, nrow(df))
  #   result[[n]] <- df
  #
  # }
  # result <- do.call(rbind, result)
  # rownames(result) <- NULL
  # return(result)



}
