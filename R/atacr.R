# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

# stop devtools::check() complain about elements in ggplot and dplyr packages
if (getRversion() >= "2.15.1")
  utils::globalVariables(c("."))

#' @importFrom magrittr %>%
#' @importFrom graphics hist
#' @importFrom stats cor kmeans median p.adjust quantile rnbinom rnorm rpois runif sd start t.test window cor.test
#' @importFrom utils capture.output read.csv str
#' @importFrom methods as
#' @importFrom SummarizedExperiment rbind
#' @importFrom stats rlnorm
no_func <-
  function(x) {
    return(FALSE)
  } #only here to make line above work

#' Get a summary of reads hitting the bait and non bait windows
#' @export
#' @param data a list of SummarizedExperiment objects from atacr::make_counts()
#' @return a table of on target and off target read counts
target_count_summary <- function(data) {
  df <- target_count_coverage(data)
  df$means <- NULL
  on_target <- off_target <- NULL #deal with devtools::check()
  d <-
    df %>% reshape::cast(sample ~ target, value = "count_sum") %>% dplyr::mutate("percent_on_target" = ((on_target /
        (
          on_target + off_target
        )) * 100))
  return(d)
}

#' Get a summary of depth of coverage in the bait and non bait windows
#' @export
#' @param data a list of SummarizedExperiment objects from atacr::make_counts()
#' @return a table of on target and off target mean depths
coverage_count_summary <- function(data) {
  df <- target_count_coverage(data)
  df$count_sum <- NULL
  return(reshape::cast(df, sample ~ target, value = "mean_coverage"))
}

#' Read count and mean coverage hitting the bait and non bait windows
#' @export
#' @param data a list of SummarizedExperiment objects from atacr::make_counts()
#' @return a dataframe of on target and off target read counts
target_count_coverage <- function(data) {
  on <- SummarizedExperiment::assay(data$bait_windows)
  off <- SummarizedExperiment::assay(data$non_bait_windows)
  target <-
    factor(c(rep("on_target", length(colnames(
      on
    ))), rep("off_target", length(colnames(
      off
    )))))
  sums <- c(colSums(on), colSums(off))
  means <- c(colMeans(on), colMeans(off))
  df <-
    data.frame(
      sample = names(sums),
      target = target,
      count_sum = sums,
      mean_coverage = means
    ) #probably not the same size?
  return(df)
}

#' identify kmeans clusters for samples
#' @export
#' @param data  a list of SummarizedExperiment objects from atacr::make_counts()
#' @param which the subdivision of the genome to calculate correlations either 'whole_genome', 'bait_windows' or 'non_bait_windows'
#' @return dataframe of cluster_id and sample name
sample_kmeans_cluster <- function(data, which = "bait_windows") {
  counts <- SummarizedExperiment::assay(data[[which]])
  k <- length(unique(data$treatments))
  c <- kmeans(t(counts), k)
  d <- data.frame(cluster_id = c$cluster)
  d$sample <- rownames(d)
  cluster_id <- NULL
  return(dplyr::arrange(d, cluster_id, sample))

}

#' count windows that have read counts below the threshold
#' @export
#' @param data  a list of SummarizedExperiment objects from atacr::make_counts()
#' @param which the subdivision of the genome to calculate correlations either
#'   'whole_genome', 'bait_windows' or 'non_bait_windows'
#' @param threshold counts windows with read counts lower than this level
#' @return dataframe of sample name, count and threshold
count_windows_under_threshold <-
  function(data,
    which = "bait_windows",
    threshold = 0) {
    counts <- SummarizedExperiment::assay(data[[which]])
    b <- apply(counts, MARGIN = 2, function(x) {
      sum(x <= threshold)
    })
    r <-
      data.frame(
        sample = names(b),
        count = b,
        threshold = rep(threshold, length(b))
      )
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
calc_quantiles <-
  function(data,
    quantiles = c(.01, .05, 0.95, 0.99),
    which = NULL) {
    if (is.null(which)) {
      bait_windows <- as.matrix(data)
      non_bait_windows <- as.matrix(data, which = "non_bait_windows")
      bwq <-
        apply(bait_windows,
          MARGIN = 2,
          quantile,
          probs = quantiles)
      non_bwq <-
        apply(non_bait_windows,
          MARGIN = 2,
          quantile,
          probs = quantiles)
      return(list(bait_windows = bwq, non_bait_windows = non_bwq))
    }
    else{
      windows <- as.matrix(data, which = which)
      return(apply(windows, MARGIN = 2, quantile, probs = quantiles))
    }
  }



se_contains_only_integers <- function(data, which) {
  a <- SummarizedExperiment::assay(data[[which]])
  return(all(a == as.integer(a)))
}



#' given a vector of values return a set of random numbers from a given
#' distribution
#' @export
#' @param obs vector of observed values
#' @param dist the distribution from which to return expected values
#' @return a vector of length obs with random variates from distribution dist
get_expected_values <- function(obs, dist = "norm") {
  exp <- rnorm(length(obs), mean = mean(obs), sd = sd(obs))
  if (dist == "pois")
    exp <- rpois(length(obs), lambda = mean(obs))
  if (dist == "nbinom") {
    est <- fitdistrplus::fitdist(obs, "nbinom")
    exp <- rnbinom(length(obs),
      size = est$estimate[['size']],
      mu = est$estimate[['mu']])
  }
  return(exp)
}

#' given a vector of numbersd returns the counts in bins of bin_width, and the count
#' @export
#' @param obs a vector of numbers
#' @param dist a string naming distribution from which to take expected counts
#' @param bin_width the width of the bins for the counts
#' @return list with members observed and expected which are vectors of counts
observed_expected_bins <-
  function(obs,
    dist = "pois",
    bin_width = 10) {
    exp <- get_expected_values(obs, dist)

    mx <- max(c(obs, exp))
    mn <- min(c(obs, exp))
    b <- seq(mn, mx + bin_width, by = bin_width)
    obs <- hist(obs, breaks = b, plot = FALSE)
    exp <- hist(exp, breaks = b, plot = FALSE)
    return(list(observed = obs$counts, expected = exp$counts))

  }

#' a median of window values across all samples in a vector, for ma plots
#' @export
#' @param sample_matrix counts extracted from a SummarizedExperiment object
#' @return the median of the provided counts, columnwise
median_virtual_experiment <- function(sample_matrix) {
  return(apply(sample_matrix, 1, median))
}

emm <- function(test, control) {
  return(log2(test) - log2(control))
}

ay <- function(test, control) {
  return(0.5 * (log2(test) + log2(control)))
}

#' given a dataframe from the estimate_fdr_multiclass() function, will return a
#' list in the format suitable for UpSetR visualisation.
#' Does not do any filtering of lists, so selected genes must be filtered before hand e.g with dplyr
#' @export
#' @param df dataframe from estimate_fdr_multiclass
#' @return list of named vectors suitable for UpSetR fromList() function
make_UpSetR <- function(df) {
  log2_fc <- direction <- a <- NULL
  r <- df %>%
    dplyr::mutate(
      direction = ifelse(log2_fc > 0, "up", "down"),
      category = paste0(direction, "_", a)
    )
  r <- r %>% split(r$category) %>%
    lapply(function(x)
      as.vector(dplyr::select(x, window)$window))
  return(r)
}

#' sim_counts - simulated count data
#'
#' The data `sim_counts` is a simulated data set with computer generated window counts for three replicates of each of two conditions in experiments with 500 bait and non-bait windows. We'll set each experiment to have 10 \% of windows differentially accessible at a difference of approximately 2 fold.
#'
#'
#' Counts in bait windows for "control" samples  will be modelled as \eqn{C \sim NB(\mu = 30, size = 10\mu)}.
#'
#' Counts in bait windows for "treatment" samples will be modelled as \eqn{C \cdot unif(0.8,1.2)}.
#'
#' Differentially accessible bait windows will be modelled as \eqn{C_{1..50} \cdot \mathcal{N}( \mu=2,\sigma = \mu/2)}
#' @format A SummarizedExperiment object
"sim_counts"

#' Simulated count data
#'
#' The data `sim_counts` is a simulated data set with computer generated window counts for three replicates of each of two conditions in experiments with 500 bait and non-bait windows. We'll set each experiment to have 10 \% of windows differentially accessible at a difference of approximately 2 fold.
#'
#' Counts in bait windows for "control" samples  will be modelled as \eqn{C \sim NB(\mu = 30, size = 10\mu)}.
#'
#' Counts in bait windows for "treatment" samples will be modelled as \eqn{C \cdot unif(0.8,1.2)}.
#'
#' Differentially accessible bait windows will be modelled as \eqn{C_{1..50} \cdot \mathcal{N}( \mu=2,\sigma = \mu/2)}
#' @format A list of SummarizedExperiment objects
"sim_counts"

#' small_counts - simulated count data
#' The data `small_counts` is basically the same thing as `sim_counts` with smaller sample of 100 bait / non-bait windows.
#' @format a list of SummarizedExperiment objects
"small_counts"

#' athal_wt_counts - real capture RNASeq count data
#' The data `athal_wt_counts` are real, experimentally derived counts from untreated WT Arabidopsis leaves for 52 baits, each set of baits representing a gene. Three replicates are provided for each gene. This data set is intended to be used in resampling procedures for making test data sets.
#' @format a named vector of counts
"athal_wt_counts"
