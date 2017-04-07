#' gets t-statistic for two vectors of data, x and y
#' @param data matrix of sample data
#' @param indices indices selected by boot::boot
#' return t the t statistic from Student's t-test
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
#'
bootstrap_t <- function(data, iterations=10){
  boot_res <- boot::boot(data, statistic=get_t, R=iterations)
  original <- boot_res$t0
  bootstraps <- boot_res$t
  p <- 1 - (sum(bootstraps > original) / iterations)
  #names(p) <- "p_val"
  return(c(original, p))
}





select_comparisons <- function(data, which = "bait_windows" , treatment_a, treatment_b){

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
  log2_fc <- apply(result, 1, function(x){ log2(x[4]) - log2(x[5])})
  return(cbind(result,log2_fc))
}



#' Estimate FDR and significantly different windows
#' @export
#'
estimate_fdr <- function(data, treatment_a, treatment_b, which = "bait_windows", iterations=10,fdr_level=0.05){

  sample_matrix <- SummarizedExperiment::assay(data[[which]])
  comparison_list <- select_comparisons(data, which, treatment_a, treatment_b)
  comparison_matrix <- cbind(comparison_list$treatment_a, comparison_list$treatment_b )
  treatment_a_names <- data$sample_names[which(data$treatments == treatment_a)]
  treatment_b_names <- data$sample_names[which(data$treatments == treatment_b)]

  #if( length(treatment_a_names) < 3 | length(treatment_b_names) < 3  ){
  #  message <- paste("Need at least 3 replicates to perform estimate bootstrap t value. Have", length(treatment_a_names), "for", treatment_a, "and", length(treatment_b_names), "for", treatment_b)
  #  stop(message)
  #}

  #if(length(treatment_a_names) != length(treatment_b_names) ){
  #  message <- paste("Must have equal number of replicates in each treatment. Have", length(treatment_a_names), "for", treatment_a, "and", length(treatment_b_names), "for", treatment_b)
  #  stop(message)
  #}

  #calc bootstrap p-values

  result <- apply(comparison_matrix, 1, bootstrap_t, iterations=iterations)
  result <- t(result)
  colnames(result) <- c("t", "p_value")
  rownames(result) <- rownames(sample_matrix)
  #calc FDR
  fdr <- p.adjust(result[,"p_value"], method="fdr")
  names(fdr) <- ("fdr")
  result <- cbind(result, fdr)

  #add means
  result <- cbind(result, get_means(comparison_list))
  #add sd
  result <- cbind(result, get_sd(comparison_list))
  #add log2 fc
  result <-  get_fc(result)


  is_sig <- result[,"fdr"] <= fdr_level


  result <- as.data.frame(result)
  result$is_sig <- is_sig
  result$window <- rownames(sample_matrix)
  #sort by fdr
  result <- dplyr::arrange(result, fdr)
  return(result[,c(10,1:9)])

}


estimate_fdr_multiclass <- function(data, common_control, which = "bait_windows", iterations=10,fdr_level=0.05){
  treatments <- data$treatments[data$treatments != common_control]
  control <- rep(common_control, length(treatments))
  comparisons <- cbind(treatments, control)
  
  r <- list()
  for (i in 1:nrow(comparisons)){
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



