#' Plot heatmap of sample count correlations
#' @export
#' @param data a list of SummarizedExperiment objects from atacr::make_counts()
#' @param which the subdivision of the genome to calculate correlations either 'whole_genome', 'bait_windows' or 'non_bait_windows'
#' @param method the correlation method to use. Any supported by `cor()` is useable
#' @return a ggplot object from geom_raster()
sample_correlation_heatmap <- function(data, which="bait_windows", method="pearson"){ # nocov start
  mat <- SummarizedExperiment::assay(data[[which]])
  hm <- get_heatmap(mat, method)
  return(hm)
}
#' generate heatmap from matrix of counts
#' @param counts a matrix of counts
#' @param method the correlation method to use, any supported by `cor()` is useable
#' @return ggplot2 plot
get_heatmap <- function(counts, method="pearson", namesort=TRUE){

  cors <- numeric(0)
  for (c1 in colnames(counts)){
    for (c2 in colnames(counts)){
      cors <- c(cors, cor(counts[,c1], counts[,c2], method=method))
    }
  }
  sample_1 <- NULL
  sample_2 <- NULL #deal with devtools::check()
  df <- expand.grid(sample_1=colnames(counts), sample_2=colnames(counts))
  df <- data.frame(sample_1=df$sample_1, sample_2=df$sample_2, correlation=cors)

  if (namesort){
   df$sample_1 <- as.character(data$sample_1)
   df$sample_1 <- factor(df$sample_1, levels=df$sample_1[order(unique(df$sample_1))])
   df$sample_2 <- as.character(data$sample_2)
   df$sample_2 <- factor(df$sample_2, levels=df$sample_2[order(unique(df$sample_2))])

  }

  correlation <- NULL #deal with devtools::check()
  heatmap <- ggplot2::ggplot(df, ggplot2::aes(sample_1, sample_2, fill = correlation)) + ggplot2::geom_raster()  + ggthemes::scale_color_ptol() + ggplot2::theme_minimal() + ggplot2::ggtitle("Sample Similarities") + ggplot2::labs(x= NULL, y=NULL)



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
  threshold <- count <- NULL #devtools::check() fix
  p <- ggplot2::ggplot(df) + ggplot2::aes(threshold, count) + ggplot2::geom_point() + ggplot2::facet_wrap(~ sample) + ggthemes::scale_color_ptol() + ggthemes::scale_fill_ptol() + ggplot2::theme_minimal()  + ggplot2::ggtitle("Counts of windows below critical threshold") + ggplot2::labs(x="Coverage threshold", y="Windows below threshold")

  return(p)
}

#' plot M (log2 ratio of a windows sample count to windows all-sample median count ) versus A (log2 sum of a windows sample count to a windows all-sample median count ) for each window
#' @export
#' @param data an atacr object
#' @param which the subset of windows to operate on
#' @param by a vector of  subset of the genome to work
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
  a <- m <- NULL
  plot <- ggplot2::ggplot(ma_df) + ggplot2::aes(a,m) + ggplot2::geom_jitter(alpha=1/length(ma_df)) + ggplot2::facet_wrap(~ sample) + ggthemes::scale_color_ptol() + ggthemes::scale_fill_ptol() + ggplot2::theme_minimal()
  return(plot)
}
#' converts SummarizedExperiment::assay matrix to a dataframe with cols 'window', 'sample' and 'count
#' @param matrix a SummarizedExperiment::assay matrix
assay_matrix_to_df <- function(matrix){
  v  <- reshape::melt( matrix )
  colnames(v) <- c("window", "sample", "count")
  return(v)
}

#' adds an 'm' and an 'a' column to an assay matrix dataframe for ma plots
#' @export
#' @param sample_matrix a SummarizedExperiment::assay from which to make the MA plot
ma_data <- function(sample_matrix){
  count <- NULL
  mve <- median_virtual_experiment(sample_matrix)
  v <- assay_matrix_to_df( sample_matrix )
  v$mve <- rep(mve, length(colnames(sample_matrix)))
  v <- dplyr::mutate(v, m = emm(count, mve))
  v <- dplyr::mutate(v, a = ay(count, mve))
  return(v)
}
#' plot the counts split by chromosome and sample
#' @param data  a list of SummarizedExperiment objects from atacr::make_counts()
#' @param which the subdivision of the genome to calculate correlations either 'whole_genome', 'bait_windows' or 'non_bait_windows'
#' @param method (bar | smooth | point) which sort of plot to return
#' @export
#' @return ggplot2 plot
plot_count_by_chromosome <- function(data, which="bait_windows", method = "bar"){
  v <- assay_matrix_to_df(SummarizedExperiment::assay(data[[which]]))
  v$window <- as.character(v$window)
  v <- tidyr::separate(v, window, into = c("seqname", "start", "stop"), sep = "[:-]")
  v$seqname <- as.factor(v$seqname)
  v$start <- as.numeric(v$start)
  v$stop <- as.numeric(v$stop)
  p <- ggplot2::ggplot(v)
  if (method == 'bar'){

    seqname <- count <- NULL
    p <- p + ggplot2::aes(start, count) + ggplot2::geom_bar(ggplot2::aes(colour=seqname, fill=seqname),stat="identity")
  }
  if (method == 'smooth'){
    p <- p + ggplot2::aes(start) + ggplot2::geom_density(ggplot2::aes(colour=seqname, fill = seqname) )
  }
  if (method == 'point'){
    p <- p + ggplot2::aes(start, count) + ggplot2::geom_point(ggplot2::aes(colour=seqname, fill = seqname) )
  }
  p <- p + ggplot2::facet_grid(sample ~ seqname) + ggthemes::scale_color_ptol() + ggthemes::scale_fill_ptol() + ggplot2::theme_minimal() + ggplot2::ggtitle("Read Count Over Chromosome") + ggplot2::labs(x= "bp", y="Read Count")
  return(p)
}

#' Plot histograms of read counts by sample and window type
#' @export
#' @param data a list of SummarizedExperiment objects from atacr::make_counts()
#' @param which the subdivision of the genome to plot (default = bait and non_bait windows)
#' @param sample the sample to plot (default = all )
#' @return a ggplot2 object
coverage_summary <- function(data, which = NULL, sample = NULL ){
  all <- as.data.frame(data)
  samp <- sample
   cov_joy_plot <- function(data){
    p <- ggplot2::ggplot(data) +
      ggplot2::aes(x = count, y =sample) +
      ggjoy::geom_joy(
        ggplot2::aes(
          fill = window_type )
      ) +
      ggplot2::facet_grid(. ~ window_type ,
        scales = "free") +
      ggthemes::scale_color_ptol()  +
      ggthemes::scale_fill_ptol() +
      ggplot2::theme_minimal() +
      ggplot2::ggtitle("Frequency of windows with read counts") +
      ggplot2::labs(x ="Read Count", y = "Number of Windows") +
      ggplot2::theme(legend.position="none")
    return(p)
  }
  if (is.null(which) & is.null(samp)){
    count <- window_type <- NULL
    p <- cov_joy_plot(all)
    return(p)
  }

  if (is.null(which) & !is.null(samp)){
    d <- all %>% dplyr::filter(sample == samp)
    p <- cov_joy_plot(d)
    return(p)
  }

  if (!is.null(which) & is.null(samp)){
    d <- all %>% dplyr::filter(window_type == which)
    p <- cov_joy_plot(d)
  return(p)
  }

  if (!is.null(which) && !is.null(samp)){
    d <- all %>% dplyr::filter(window_type == which) %>%
      dplyr::filter(sample == samp)
    p <- cov_joy_plot(d)
    return(p)
  }

}

#' Plot density of read counts by sample over the chromosomes
#' @export
#' @param data a list of SummarizedExperiment objects from atacr::make_counts()
#' @param which the subdivision of the genome to plot (default = bait and non_bait windows)
#' @return a ggplot2 object
chromosome_coverage <- function(data, which = NULL){
  all <- as.data.frame(data)
  d <- NULL
  if (is.null( which )){
    d <- all
  }
  else{
    window_type <- NULL #deal with devtools::check()
    d <- all %>% dplyr::filter(window_type == which)
  }
  p <- ggplot2::ggplot(d) + ggplot2::aes(start) + ggplot2::geom_density(ggplot2::aes(colour=window_type)) + ggplot2:: facet_grid( sample ~ chromosome, scales = "free_x" ) + ggthemes::scale_color_ptol() + ggthemes::scale_fill_ptol() + ggplot2::theme_minimal()  + ggplot2::ggtitle("Density of coverage over chromosomes")
  return(p)
}

#' Named distribution qqplot
#' @export
#' @param obs observed values
#' @param dist expected distribution
#' @return ggplot2 object
qqarb <- function(obs, dist = "norm" ){
  exp <- get_expected_values(obs, dist)

  df <- data.frame(observed_values = sort(obs), expected_values =  sort(exp) )
  expected_values <- observed_values <- NULL
  p <- ggplot2::ggplot(df) +ggplot2::aes(expected_values, observed_values) + ggplot2::geom_point() + ggplot2::geom_abline(intercept = 0, slope = 1) + ggthemes::scale_color_ptol() + ggthemes::scale_fill_ptol() + ggplot2::theme_minimal()
  return(p)
}

get_mart <- function(ensembl, ens_dataset){
  mart <- switch(ensembl,
    plants = biomaRt::useMart("plants_mart",
      host="plants.ensembl.org",
      dataset=ens_dataset),
    ensembl = biomaRt::useMart("ensembl", dataset=ens_dataset)
  )
  return(mart)
}

get_gene_coords <- function(gene_id, mart){
  filter <- switch(mart@biomart,
    plants_mart = "tair_locus",
    ENSEMBL_MART_ENSEMBL = "with_entrezgene"
    )

  coords <- biomaRt::getBM(attributes=c("chromosome_name", "start_position", "end_position", "strand"), filters=filter, values = gene_id, mart=mart )
  return(coords)
}
select_colours <-function(data){
  t <- treatments(data)
  allcols <- RColorBrewer::brewer.pal(8, "Dark2")
  allcols <- rep(allcols, ceiling(length(t) / length(allcols)) )
  cols <- allcols[1:length(t)]
  names(cols) <- t
  tr <- NULL
  for (i in data$sample_names){tr <- c(tr, treatment_from_name(data, i))}
  return(cols[tr])
}

get_coverage_in_regions <- function(data, which, coords){
  sname <- coords$chromosome_name[[1]]
  strt <- coords$start_position[[1]]
  stp <- coords$end_position[[1]]
  strand <- coords$strand[[1]]

  roi <- GenomicRanges::GRanges(seqnames=sname, ranges=strt:stp )
  se <- IRanges::subsetByOverlaps(data[[which]], roi)
  colrs <- select_colours(data) #c(rep("grey", 2), rep("red", 2))
  intens <- GenomeGraphs::makeGenericArray(
    intensity = SummarizedExperiment::assay(se),
    probeStart = se@rowRanges@ranges@start,
    #probeEnd = (se@rowRanges@ranges@start + se@rowRanges@ranges@width),
#    nProbes = nrow(se),
 #   probeId = se@rowRanges@ranges@NAMES,
    dp = GenomeGraphs::DisplayPars(color = colrs, size=2, lwd = 2, type="line", pointSize = 1)
  )
  return(intens)
}

#' coverage over gene model
#' @export
#' @param data atacr object
#' @param gene_id the id of the gene to plot around
#' @param which the subset of the data to plot.
#' @param ensembl one of 'plants', 'ensembl' - which version of ensembl to connect to
#' @param ens_dataset which ensembl dataset to connect to
#' @return plot object
view_gene <- function(data,gene_id,which="bait_windows", ensembl = "plants", ens_dataset = "athaliana_eg_gene"){

  ##get connection to biomart
  mart <- get_mart(ensembl, ens_dataset)

  ##extract gene coords from ensembl
  coords <- get_gene_coords(gene_id, mart)
  start <- coords$start_position[[1]]
  end <- coords$end_position[[1]]
  chrom <- coords$chromosome[[1]]
  strand <- as.character(coords$strand[[1]])
  ##get coverage count in each window over coords
  values <- get_coverage_in_regions(data, which, coords)

  if(strand ==  "1"){
    axis <- GenomeGraphs::makeGenomeAxis(add53 = TRUE, littleTicks=TRUE, dp =GenomeGraphs::DisplayPars(cex=0.5) )
    strand = "+"
  } else {
    axis <- GenomeGraphs::makeGenomeAxis(add35 = TRUE, littleTicks=TRUE, dp =GenomeGraphs::DisplayPars(cex=0.5))
    strand = "-"
    }

  ##get features in gene region from ensembl

  g <- GenomeGraphs::makeGeneRegion(
    start = start,
    end = end,
    chromosome = chrom,
    strand = strand,
    biomart = mart,
    dp = GenomeGraphs::DisplayPars(protein_coding = "steelblue")
    )

  view <- list(
    GenomeGraphs::makeTitle(gene_id),
    "counts" = values,
    "gene" = g,
    axis
    )

  GenomeGraphs::gdPlot(view, minBase = start, maxBase = end)


}
# nocov end
