#' @export
all_versus_all_normalized_plot <- function(df){
  col_count <- dplyr::n_distinct(df$Sample) + 1
  p <- df %>%
    dplyr::select(Bait, Sample, log_normalized_value) %>%
    dplyr::mutate(log_normalized_value = ifelse(is.na(log_normalized_value),0, log_normalized_value)) %>%
    dplyr::mutate(log_normalized_value = ifelse(is.infinite(log_normalized_value),0, log_normalized_value)) %>%
    tidyr::spread(Sample, log_normalized_value)  %>%
    GGally::ggpairs(columns=2:col_count, lower=list(continuous="points"), upper=list(continuous=GGally::wrap("cor", size = 10)))
  return(p)
}

#' @export
all_versus_all_read_count_plot <- function(df){
  col_count <- dplyr::n_distinct(df$Sample) + 1
  p <- df %>%
    dplyr::select(Bait, Sample, sample_bait_read_count) %>%
    dplyr::mutate(sample_bait_read_count = ifelse(is.na(sample_bait_read_count),0, sample_bait_read_count)) %>%
    dplyr::mutate(sample_bait_read_count = ifelse(is.infinite(sample_bait_read_count),0, sample_bait_read_count)) %>%
    tidyr::spread(Sample, sample_bait_read_count)  %>%
    GGally::ggpairs(columns=2:col_count, lower=list(continuous="points"), upper=list(continuous=GGally::wrap("cor", size = 10)))
  return(p)
}


#' Plot on target and off target reads
#' @export
#' @param data a list of SummarizedExperiment objects from atacr::make_counts()
#' @param scales, either "free_y" for free y axes per subplot, or "fixed" for fixed ones.
#' @return a ggplot barchart object
on_target_plot <- function(data, scales="free_y"){
  df <- on_target_count(data)
  p <- ggplot2::ggplot(df) + ggplot2::aes(target, count_sum) + ggplot2::geom_bar(ggplot2::aes(colour=target, fill=target), stat="identity", ) + ggplot2::facet_wrap( ~ sample, scales="free_y")
  return(p)
}

#' Plot heatmap of sample count correlations
#' @export
#' @param data a list of SummarizedExperiment objects from atacr::make_counts()
#' @param which the subdivision of the genome to calculate correlations either 'whole_genome', 'bait_windows' or 'non_bait_windows'
#' @return a ggplot object from geom_raster()
sample_correlation_heatmap <- function(data, which="bait_windows", method="pearson"){
  mat <- SummarizedExperiment::assay(data[[which]])
  hm <- get_heatmap(mat, method)
  return(hm)
}
#' generate heatmap from matrix of counts
#' @param mat a matrix of counts
get_heatmap <- function(counts,method="pearson"){

  cors <- numeric(0)
  for (c1 in colnames(counts)){
    for (c2 in colnames(counts)){
      cors <- c(cors, cor(counts[,c1], counts[,c2], method=method))
    }
  }

  df <- expand.grid(sample_1=colnames(counts), sample_2=colnames(counts))

  df <- data.frame(sample_1=df$sample_1, sample_2=df$sample_2, correlation=cors)
  heatmap <- ggplot2::ggplot(df, ggplot2::aes(sample_1, sample_2, fill = correlation)) + ggplot2::geom_raster()
  return(heatmap)
}
#' generate cumulative plot of number of windows below a threshold in samples
#' @export
windows_below_coverage_threshold_plot <- function(data, which="bait_windows", from=0, to=10){
  df <- count_windows_under_threshold(data, which = which, threshold = from)
  for (i in (from + 1):to ){
    df <- rbind(df, count_windows_under_threshold(data, which = which, threshold = i))
  }
  rownames(df) <- NULL

  a <- df %>% dplyr::group_by(sample) %>% dplyr::mutate(windows_below_threshold = cumsum(count)) %>% dplyr::arrange(sample, threshold)
  p <- ggplot2::ggplot(a) + ggplot2::aes(threshold, windows_below_threshold) + geom_point() + ggplot2::facet_wrap(~ sample)

  return(p)
}

#' @export
ma_plot <- function(data, which = "bait_windows", by = "none" ){
  mve <- median_virtual_experiment(data, which = which)
  v  <- apply(SummarizedExperiment::assay(data[[which]]), 2, ma, control = mve)
  return(v) ## doesnt work!! need to draw ma plots. Perhaps avoid using apply to avoid lists...
}
