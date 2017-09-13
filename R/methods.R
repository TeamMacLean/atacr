

meta_summary <- function(atcr) {
  samples = paste(unique(atcr$sample_names), collapse = ",")
  treatments = paste(unique(atcr$treatments), collapse = ",")
  sample_count = length(unique(atcr$sample_names))
  treat_count  = length(unique(atcr$treatments))
  return(
    cat(
      "ATAC-seq experiment of",
      treat_count,
      "treatments in",
      sample_count,
      "samples\n",
      "Treatments:",
      treatments,
      "\n",
      "Samples:",
      samples,
      "\n",
      "Bait regions used:",
      length(atcr$bait_regions),
      "\n",
      "Total Windows:",
      length(atcr$whole_genome) ,
      "\n"

    )
  )
}

#' writes a summary of the metadata for a given atacr object
#' @export
#' @param x an atacr object
#' @param \dots other options for print generic
print.atacr <- function(x, ...) {
  meta_summary(x)
  invisible(x)
}

#' writes a detailed data summary of the atacr object
#' @export
#' @param object an atacr object
#' @param \dots other options for summary generic
summary.atacr <- function(object, ...) {
  atcr <- object
  meta <- meta_summary(atcr)
  on_target <-
    paste(capture.output(target_count_summary(atcr)), collapse = "\n")
  coverage <-
    paste(capture.output(coverage_count_summary(atcr)), collapse = "\n")
  quantiles <-
    paste(capture.output(calc_quantiles(atcr)), collapse = "\n")
  return(
    cat(
      meta,
      "\n",
      "On/Off target read counts:\n",
      on_target,
      "\n",
      "Quantiles:",
      "\n",
      quantiles,
      "\n",
      "Read depths:\n",
      coverage
    )
  )

}
#' returns given subset of data in atacr object as a matrix
#' @export
#' @param x an atacr object
#' @param \dots other options for generic
#' @param which the subset of data to work on
#' @return matrix of counts in subset
as.matrix.atacr <- function(x, ..., which = "bait_windows") {
  atcr <- x
  return(SummarizedExperiment::assay(atcr[[which]]))
}

#' returns dataframe of data in atacr object
#' @export
#' @param x object to print
#' @param \dots other options for generic
#' @return dataframe
as.data.frame.atacr <- function(x, ...) {
  atcr <- x
  if (is.null(atcr[["dataframe"]])) {
    bw <- as.matrix.atacr(atcr, which = "bait_windows")
    nbw <- as.matrix.atacr(atcr, which = "non_bait_windows")
    bw_df <- reshape2::melt(bw)
    colnames(bw_df) <- c("name", "sample", "count")
    bw_df$window_type <- factor(rep("bait_windows", nrow(bw_df)))
    nbw_df <- reshape2::melt(nbw)
    colnames(nbw_df) <- c("name", "sample", "count")
    nbw_df$window_type <-
      factor(rep("non_bait_windows", nrow(nbw_df)))
    df <- rbind(bw_df, nbw_df)
    df$name <- stringr::str_replace(df$name, "-$", "minus")
    df$name <- stringr::str_replace(df$name, "\\+$", "plus")
    name <- NULL #deal with NSE of devtools::check()
    df <-
      tidyr::separate(df, name, c("chromosome", "start", "stop", "strand"), sep =
          '[-:]')
    df$start <- as.integer(df$start)
    df$stop <- as.integer(df$stop)
    df$chromosome <- factor(df$chromosome)
    atcr[["dataframe"]] <- df
    return(df)
  }
  else{
    return(atcr[["dataframe"]])
  }

}

#' returns summary plot of data in atacr object
#' @method plot atacr
#' @S3method plot atacr
#' @param x atacr object
#' @param \dots extra options for generic
#' @return gridExtra plot
plot.atacr <- function(x, ...) {
  atcr <- x
  #histogram of coverages by sample and window type
  p1 <- coverage_summary(atcr)

  #density of coverage by chromosome region, bait windows
  p2 <-  chromosome_coverage(atcr)

  return(gridExtra::grid.arrange(p1, p2, nrow = 2))

}
