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
      print(w)
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
#' @return vector of 2 items, observed value t statisitc and p, calculated as proportion of bootstrap iterations greater than original t
bootstrap_t <- function(data, iterations=10){
  boot_res <- boot::boot(data, statistic=get_t, R=iterations)
  original <- boot_res$t0
  bootstraps <- boot_res$t
  p <- (sum(bootstraps > original) / iterations)
  if ( is.nan(original) | is.na(original)){
    p <- original
  }

  else if (original < 0) {
    p <- sum(bootstraps < original) / iterations
  }
  #names(p) <- "p_val"
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

get_means <- function(comparisons){

  mean_count_a <- apply(comparisons$treatment_a_data, 1, mean)
  mean_count_b <- apply(comparisons$treatment_b_data, 1, mean)
  return(cbind(mean_count_a,mean_count_b))
}

get_sd <- function(comparisons){
  sd_a <- apply(comparisons$treatment_a_data, 1, sd)
  sd_b <- apply(comparisons$treatment_b_data, 1, sd)
  return(cbind(sd_a,sd_b))
}

get_fc <- function(result){
  mean_count_a <- mean_count_b <- NULL
  result <- dplyr::mutate(result, log2_fc = log2(mean_count_a) - log2(mean_count_b))
  return(result)
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
      treatment_b_names = data$sample_names[which(data$treatments == treatment_b)]
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
#' @param remove_zeroes remove all windows that have a summed count of zero before computing
estimate_fdr <- function(data, treatment_a, treatment_b, which = "bait_windows", iterations=10,fdr_level=0.05, remove_zeroes = FALSE){


  d <- select_data(data, treatment_a, treatment_b, which)
  check_data(d, treatment_a, treatment_b)

  if (remove_zeroes) {
    x <- rowSums(d$counts ) > 0
    d$counts <- d$counts[x, ]

  }
  #calc bootstrap p-values

  result <- apply(d$counts, 1, bootstrap_t, iterations = iterations)
  result <- t(result)
  colnames(result) <- c("t", "p_value")

  result <- data.frame(
    window = rownames(result),
    t = result[, "t"],
    p_value = result[,"p_value"]
  )

  #rownames(result) <- rownames(d$counts)
  #calc FDR
  result$fdr <- p.adjust(result[,"p_value"], method = "fdr")
  #names(fdr) <- ("fdr")
  #result <- cbind(result, fdr)

  is_sig <- result[,"fdr"] <= fdr_level


  result <- as.data.frame(result)
  result$is_sig <- is_sig
  result$window <- rownames(d$counts)

  #add means
  result <- cbind(result, get_means(d$comparisons))
  #add sd
  result <- cbind(result, get_sd(d$comparisons))
  #add log2 fc
  result <-  get_fc(result)

  #sort by fdr
  #result <- dplyr::arrange(result, fdr)
  return(result)#[,c(10,1:9)])

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
#' @param remove_zeroes remove all windows that have a summed count of zero before computing
#' @return a dataframe
estimate_bayes_factor <- function(atacr, treatment_a, treatment_b, which = "bait_windows", factor = 4, remove_zeroes= FALSE){

  d <- select_data(atacr, treatment_a, treatment_b, which)
  check_data(d, treatment_a, treatment_b)

  if (remove_zeroes) {
    x <- rowSums(d$counts ) > 0
    d$counts <- d$counts[x, ]
  }


  result <- apply(as.data.frame(d$counts), 1, bayes_t, treatment_a_names = d$treatment_a_names, treatment_b_names = d$treatment_b_names)


  result <- data.frame(
    window = names(result),
    bayes_factor = result
    )
  result$is_sig <- (result$bayes_factor >= factor)
  result <- cbind(result, get_means(d$comparisons))
  #add sd
  result <- cbind(result, get_sd(d$comparisons))
  #add log2 fc
  result <-  get_fc(result)


  return(result)
}

#' Estimate BayesFactor and mark significantly different windows for many experiments
#' @export
#' @param data an atacr object
#' @param common_control the treatment to consider the control for all other treatments
#' @param which the subset of windows to consider
#' @param factor the BayesFactor to consider significant
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
    r[[i]] <- df
  }
  return(do.call(rbind, r))

}
#' Estimate differential window counts  and mark significantly different windows using edgeR exact method for two samples
#' @export
#' @param data an atacr object
#' @param common_control the treatment to consider the control for all other treatments
#' @param which the subset of windows to consider
#' @param sig_level the p_value to consider significant
edgeR_exact <- function(data, which = "bait_windows", treatment_a = NULL, treatment_b = NULL, remove_zeros = FALSE, sig_level = 0.05 ){

  dlist <- select_data(data, treatment_a, treatment_b, which)

  group <- c(rep(treatment_a, length(dlist$treatment_a_names)), rep(treatment_b, length(dlist$treatment_b_names)) )

  dg <- edgeR::DGEList(dlist$counts, group = group, remove.zeros = remove_zeros)
  dg <- edgeR::estimateDisp(dg)
  et <- edgeR::exactTest(dg)
  names <- rownames(et$table)

  result <- data.frame(
    window = rownames(et$table),
    p_value = et$table$PValue
  )
  result$is_sig <- (result$p_value <= sig_level)
  result <- cbind(result, get_means(dlist$comparisons))
  #add sd
  result <- cbind(result, get_sd(dlist$comparisons))
  #add log2 fc
  result <-  get_fc(result)
  return(result)
}
#' Estimate differential window counts  and mark significantly different windows using edgeR glmFIT method for multiple samples with common control
#' @export
#' @param data an atacr object
#' @param treatment_a the first treatment to consider
#' @param treatment_b the second treatment to consider
#' @param which the subset of windows to consider
#' @param sig_level the p_value to consider significant
edgeR_multiclass <- function(data, common_control, which = "bait_windows", sig_level = 0.05, remove.zeros = FALSE){

  ctrl_idcs <- which(data$treatments == common_control)
  other_idcs <- which(data$treatments != common_control)
  new_order <- c(ctrl_idcs, other_idcs)
  treatments <- as.factor(data$treatments[ new_order ])
  samples <- data$sample_names[ new_order ]
  df <- data.frame(sample = samples, treatment = as.factor(as.numeric(treatments)))
  design <- model.matrix(~treatment, data = df)
  num_levels <- nlevels(as.factor(treatments))
  dglist <- edgeR::DGEList(SummarizedExperiment::assay(data[[which]]), remove.zeros = remove.zeros)

  dglist <- edgeR::estimateDisp(dglist, design)
  fit <- edgeR::glmQLFit(dglist, design)

  dgelrts <- list()

  for (i in 2:num_levels) {
    curr_t <- levels(treatments)[i]
    dgelrts[[curr_t]] <- edgeR::glmQLFTest(fit, coef = i)
  }

  result <- list()

  for(n in names(dgelrts)){
    tb <- dgelrts[[n]]$table
    df <-  data.frame(
      window = rownames(tb),
      p_value = tb$PValue,
      f = tb$F
    )

    dlist <- select_data(data, n, common_control, which)
    df$is_sig <- (df$p_value <= sig_level)
    df <- cbind(df, get_means(dlist$comparisons))

    #add sd
    df <- cbind(df, get_sd(dlist$comparisons))
    #add log2 fc
    df <- get_fc(df)
    df$a <- rep(n, nrow(df))
    df$b <- rep(common_control, nrow(df))
    result[[n]] <- df
  }
  result <- do.call(rbind, result)
  rownames(result) <- NULL
  return(result)



}
