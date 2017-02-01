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
#' @return ggplot2 plot
get_heatmap <- function(counts){

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
#' @param data a list of SummarizedExperiment objects from atacr::make_counts()
#' @param which ("bait_windows") the subdivision of the genome to calculate correlations either 'whole_genome', 'bait_windows' or 'non_bait_windows'
#' @param from (0) the lowest threshold to consider
#' @param to (10) the highest threshold to consider
#' @export
#' @return ggplot2 plot
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

#' plot M (log2 ratio of a windows sample count to windows all-sample median count ) versus A (log2 sum of a windows sample count to a windows all-sample median count ) for each window
#' @export
ma_plot <- function(data, which = "bait_windows", by = NULL ){
  sample_matrix <- matrix(0)
 # by is to decide on sub-group, IE whole window, chromosome, region
  if (!is.null(by)){
    roi <- GenomicRanges::GRanges(seqnames = by)
    sample_matrix <- SummarizedExperiment::assay(IRanges::subsetByOverlaps(data[[which]],roi))

  }
  else{
    #print(colnames(SummarizedExperiment::assay(data[[which]])))
    #print(which)
    sample_matrix <- SummarizedExperiment::assay(data[[which]])
    #print(colnames(sample_matrix))
    #print(str(sample_matrix))
  }
  ma_df <-  ma_data(sample_matrix)
  #print(str(ma_df))
  #  do ggplot
  plot <- ggplot2::ggplot(ma_df) + ggplot2::aes(a,m) + geom_jitter() + facet_wrap(~ sample)
  return(plot)
}
#' converts SummarizedExperiment::assay matrix to a dataframe with cols 'window', 'sample' and 'count
assay_matrix_to_df <- function(matrix){
  v  <- reshape::melt( matrix )
  colnames(v) <- c("window", "sample", "count")
  return(v)
}

#' adds an 'm' and an 'a' column to an assay matrix dataframe for ma plots
ma_data <- function(sample_matrix){
  mve <- atacr::median_virtual_experiment(sample_matrix)
  v <- atacr::assay_matrix_to_df( sample_matrix )
  v$mve <- rep(mve, length(colnames(sample_matrix)))
  v <- dplyr::mutate(v, m = atacr::emm(count, mve))
  v <- dplyr::mutate(v, a = atacr::ay(count, mve))
  return(v)
}
#' plot the counts split by chromosome and sample
#' @param data  a list of SummarizedExperiment objects from atacr::make_counts()
#' @param which the subdivision of the genome to calculate correlations either 'whole_genome', 'bait_windows' or 'non_bait_windows'
#' @param method (bar | smooth | point) which sort of plot to return
#' @export
#' @return ggplot2 plot
plot_count_by_chromosome <- function(data, which="bait_windows", method = "bar"){
  v <- atacr::assay_matrix_to_df(SummarizedExperiment::assay(data[[which]]))
  v$window <- as.character(v$window)
  v <- tidyr::separate(v, window, into = c("seqname", "start", "stop"), sep = "[:-]")
  v$seqname <- as.factor(v$seqname)
  v$start <- as.numeric(v$start)
  v$stop <- as.numeric(v$stop)
  p <- ggplot2::ggplot(v)
  if (method == 'bar'){

    p <- p + ggplot2::aes(start, count) + ggplot2::geom_bar(ggplot2::aes(colour=seqname, fill=seqname),stat="identity")
  }
  if (method == 'smooth'){
    p <- p + ggplot2::aes(start) + ggplot2::geom_density(ggplot2::aes(colour=seqname, fill = seqname) )
  }
  if (method == 'point'){
    p <- p + ggplot2::aes(start, count) + ggplot2::geom_point(ggplot2::aes(colour=seqname, fill = seqname) )
  }
  p <- p + facet_grid(sample ~ seqname)
  return(p)
}
