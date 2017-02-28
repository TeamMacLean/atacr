# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

#' @importFrom magrittr %>%
no_func <- function(x){return(FALSE)} #only here to make line above work


#' Get a summary of reads hitting the bait and non bait windows
#' @export
#' @param data a list of SummarizedExperiment objects from atacr::make_counts()
#' @return a table of on target and off target read counts
target_count_summary <- function(data){
  df <- target_count_coverage(data)
  df$means <- NULL
  return(reshape::cast(df, sample ~ target, value="count_sum"))
}
#' Get a summary of depth of coverage in the bait and non bait windows
#' @export
#' @param data a list of SummarizedExperiment objects from atacr::make_counts()
#' @return a table of on target and off target mean depths
coverage_count_summary <- function(data){
  df <- target_count_coverage(data)
  df$count_sum <- NULL
  return(reshape::cast(df, sample ~ target, value="mean_coverage"))
}


#' Read count and mean coverage hitting the bait and non bait windows
#' @export
#' @param data a list of SummarizedExperiment objects from atacr::make_counts()
#' @return a dataframe of on target and off target read counts
target_count_coverage <- function(data){
  on <- SummarizedExperiment::assay(data$bait_windows)
  off <- SummarizedExperiment::assay(data$non_bait_windows)
  target <- factor( c( rep("on_target", length(colnames(on))), rep("off_target", length(colnames(off)))))
  sums <- c(colSums(on), colSums(off) )
  means <- c(colMeans(on), colMeans(off))
  df <- data.frame(sample = names(sums), target = target, count_sum = sums, mean_coverage = means) #probably not the same size?
  return(df)
}



#' identify kmeans clusters for samples
#' @export
#' @param data  a list of SummarizedExperiment objects from atacr::make_counts()
#' @param which the subdivision of the genome to calculate correlations either 'whole_genome', 'bait_windows' or 'non_bait_windows'
#' @return dataframe of cluster_id and sample name
sample_kmeans_cluster <- function(data, which="bait_windows"){
  counts <- SummarizedExperiment::assay(data[[which]])
  k <- length(unique(data$treatments))
  c <- kmeans(t(counts), k)
  d <- data.frame(cluster_id = c$cluster)
  d$sample <- rownames(d)
  return(dplyr::arrange(d, cluster_id, sample))

}

#' count windows that have read counts below the threshold
#' @export
#' @param data  a list of SummarizedExperiment objects from atacr::make_counts()
#' @param which the subdivision of the genome to calculate correlations either 'whole_genome', 'bait_windows' or 'non_bait_windows'
#' @param threshold counts windows with read counts lower than this level
#' @return dataframe of sample name, count and threshold
count_windows_under_threshold <- function(data, which="bait_windows", threshold=0){
  counts <- SummarizedExperiment::assay(data[[which]])
  b <- apply(counts, MARGIN = 2, function(x){sum(x <= threshold)})
  r <- data.frame(sample=names(b), count=b, threshold=rep(threshold, length(b)) )
  rownames(r) <- NULL
  return(r)
}


#' a median of window values across all samples in a vector, for ma plots
#' @param data a list of SummarizedExperiment objects from atacr::make_counts()
median_virtual_experiment <- function(sample_matrix ){
  return(apply(sample_matrix, 1, median))
}

emm <- function(test,control){
  return(log2(test) - log2(control) )
}

ay <- function(test,control){
  return(0.5 * (log2(test) + log2(control)))
}
