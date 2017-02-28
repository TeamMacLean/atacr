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
  log2_fold_change <- apply(result, 1, function(x){ log2(x[4]) - log2(x[5])})
  return(cbind(result,log2_fold_change))
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
  print(head(result))
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
  #sort by fdr
  result <- dplyr::arrange(result, fdr)
  return(result)

}

estimate_rp <- function(treatment_a, treatment_b, read_exp_mapping, data, iterations=10, fdr_level=0.05){



}

## how this is supposed to work:
# 1. create file of treatment, sample_name, bam_file_path to load in with read_experiment_info
# this returns object like
# 2. df <- data.frame(treatment=c("mock", "flg22", "UB40","mock","mock"), "sample_name"=c("sample_1", "atac_3", "sample_4", "atac_2", "atac_5"), "bam_file_path"=c("result/bam.bam", "other/bam.sorted.bam", "last/bam.bam.bam", "bam", "no.bam"))
# 3. Then you can get a list of names for a sample from the df using
# n <- names_from_treatment(c("mock"), df)
# 4. Then you can subset the experiment matrix e.g in the RangedSummarizedExperiment object by name (nb you'll need to have set colnames() to the sample_names in the df)
# mock_experiments <- assay(filtered)[,n]


