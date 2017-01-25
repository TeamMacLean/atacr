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
