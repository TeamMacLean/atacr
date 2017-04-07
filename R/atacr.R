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
  d <- df %>% reshape::cast( sample ~ target, value="count_sum") %>% dplyr::mutate("percent_on_target" = ((on_target/(on_target + off_target)) * 100) )
  return(d)
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
#' @param which the subdivision of the genome to calculate correlations either
#'   'whole_genome', 'bait_windows' or 'non_bait_windows'
#' @param threshold counts windows with read counts lower than this level
#' @return dataframe of sample name, count and threshold
count_windows_under_threshold <- function(data, which="bait_windows", threshold=0){
  counts <- SummarizedExperiment::assay(data[[which]])
  b <- apply(counts, MARGIN = 2, function(x){sum(x <= threshold)})
  r <- data.frame(sample=names(b), count=b, threshold=rep(threshold, length(b)) )
  rownames(r) <- NULL
  return(r)
}

#' report counts at each quantile for each sample
#' @export
#' @param data a list of SummarizedExperiment objects from atacr::make_counts()
#' @param quantiles a vector of quantiles to report
#' @param which the subset of data windows to report on. Default =
#'   "bait_windows" and "non_bait_windows"
#' @return list of counts at quantiles
calc_quantiles <- function(data, quantiles = c(.01,.05,0.95, 0.99), which = NULL){
  if( is.null(which) ){
    bait_windows <- as.matrix(data)
    non_bait_windows <- as.matrix(data, which = "non_bait_windows")
    bwq <- apply(bait_windows, MARGIN = 2, quantile, probs = quantiles )
    non_bwq <- apply(non_bait_windows, MARGIN = 2, quantile, probs = quantiles )
    return(list(bait_windows = bwq, non_bait_windows = non_bwq))
  }
  else{
    windows <- as.matrix(data, which = which)
    return( apply(windows, MARGIN = 2, quantile, probs = quantiles ) )
  }
}

#' run fitdistrplus to estimate chisquared goodness of fit and kolmogorov
#' smirnov tests for fit to distributions
#' @export
#' @param dist a distribution name recognised by fitdistrplus
#' @param data a list of SummarizedExperiment objects from atacr::make_counts()
#' @param keep a character vector of fitdistrplus statistics to report
#' @param which the subset of data windows to report on. Default =
#'   "bait_windows"
#' @return data.frame with columns for keep statistics, dsitribution tested and
#'   samples
get_fit <- function(dist, data = NULL, which = "bait_windows", keep = c("chisqpvalue", "cvm", "ad", "ks"), min_count = 0){
  d <- as.matrix(data, which = which)
  d[d < min_count] <- NA
  r <- apply(d, 2, function(x, dist) fitdistrplus::fitdist(x[!is.na(x)], dist), dist )
  gof <- lapply(r, fitdistrplus::gofstat)

  r<- matrix(0, nrow=length(names(gof)), ncol=length(keep) )
  colnames(r) <- keep
  l <- lapply(gof, function(x) x[keep])
  rownames(r) <- names(l)
  num_rows <-length(names(l))
  for (row in 1:num_rows){
    for (col in 1:length(names(l[[row]]))){
      val <- unname(l[[row]][[col]])
      if (length(val) == 0){ val <- NA }
      r[row,col] <- val
    }
  }
  r <- as.data.frame(r)
  r$sample <- rownames(r)
  r$distribution <- rep(dist, length(r$sample))
  rownames(r) <- NULL
  return(r)
}

#' run distribution fitting for all samples in a data set using multiple named
#' distributions
#' @param data a list of SummarizedExperiment objects from atacr::make_counts()
#' @param dists a character vector of distribution names recognised by
#'   fitdistrplus
#' @param which the subset of data windows to report on. Default =
#'   "bait_windows"
#' @param keep a character vector of fitdistrplus statistics to report
#' @return data.frame with columns for keep statistics,  and samples
get_fits <- function(data, 
                     dists = c("norm", "pois","nbinom"), 
                     which = "bait_windows", 
                     keep = c("chisqpvalue", "cvm", "ad", "ks"), 
                     min_count = 0 ){
  
  result <- lapply(dists, get_fit, 
                   data, 
                   which = which, 
                   keep = keep, 
                   min_count = min_count
                   )
  return(do.call("rbind", result))
}

#' given a vector of values return a set of random numbers from a given
#' distribution
#' @export
#' @param obs vector of observed values
#' @param dist the distribution from which to return expected values
#' @return a vector of length obs with random variates from distribution dist
get_expected_values <- function(obs, dist="norm"){
  exp <- rnorm(length(obs), mean = mean(obs), sd = sd(obs))
  if (dist == "pois") exp <- rpois(length(obs), lambda = mean(obs) )
  if (dist == "nbinom"){
    est <- fitdistrplus::fitdist(obs, "nbinom")
    exp <- rnbinom(length(obs), 
                   size = est$estimate[['size']], 
                   mu = est$estimate[['mu']]
                   )
  }
  return(exp)
}

#' given a vector of numbersd returns the counts in bins of bin_width, and the count
#' @export
#' @param obs a vector of numbers
#' @param dist a string naming distribution from which to take expected counts
#' @param bin_width the width of the bins for the counts
#' @return list with members observed and expected which are vectors of counts
observed_expected_bins <- function(obs, dist = "pois", bin_width=10) {
  exp <- get_expected_values(obs,dist)

  mx <- max(c(obs,exp))
  mn <- min(c(obs,exp))
  b <- seq(mn, mx + bin_width, by=bin_width)
  obs <- hist(obs, breaks = b, plot=FALSE)
  exp <- hist(exp, breaks = b, plot=FALSE)
  return(list(observed=obs$counts, expected=exp$counts))

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

#' given a dataframe from the estimate_fdr_multiclass() function, will return a 
#' list in the format suitable for UpSetR visualisation.
#' Does not do any filtering of lists, so selected genes must be filtered before hand e.g with dplyr
#' @export
#' @param df dataframe from estimate_fdr_multiclass
#' return list of named vectors suitable for UpSetR fromList() function
make_UpSetR <- function(df){
  r <- df %>% 
    dplyr::mutate(
      direction = ifelse(log2_fc > 0, "open", "closed"), 
      category = paste0(direction, "_", a)
      ) 
  r <- r %>% split( r$category  ) %>% 
    lapply( function(x) as.vector(dplyr::select(x, window)$window))
  return(r)
}
