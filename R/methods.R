meta_summary <- function(atcr){
  samples = paste(unique(atcr$sample_names), collapse=",")
  treatments = paste(unique(atcr$treatments), collapse=",")
  sample_count = length(unique(atcr$sample_names))
  treat_count  = length(unique(atcr$treatments))
  return(cat(
    "ATAC-seq experiment of", treat_count, "treatments in", sample_count, "samples\n",
    "Treatments:", treatments, "\n",
    "Samples:", samples, "\n",
    "Bait regions used:", length(atcr$bait_regions), "\n",
    "Total Windows:", length(atcr$whole_genome) , "\n"

    ))
}

print.atacr <- function(atcr){
  meta_summary(atcr)
}

summary.atacr <- function(atcr){
  meta <- meta_summary(atcr)
  on_target <- paste(capture.output(target_count_summary(atcr)), collapse="\n")
  coverage <- paste(capture.output(coverage_count_summary(atcr)), collapse="\n")
  quantiles <- paste(capture.output(calc_quantiles(atcr)), collapse="\n" )
  return(cat(
    meta, "\n",
    "On/Off target read counts:\n",
    on_target,"\n",
    "Quantiles:", "\n",
    quantiles, "\n",
    "Read depths:\n",
    coverage
    ) )


}

as.matrix.atacr <- function(atcr, which = "bait_windows"){
  return( SummarizedExperiment::assay(atcr[[which]]))
}

as.data.frame.atacr <- function(atcr){
  if (is.null(atcr[["dataframe"]])){
    bw <- as.matrix.atacr(atcr, which = "bait_windows")
    nbw <- as.matrix.atacr(atcr, which = "non_bait_windows")
    bw_df <- reshape2::melt(bw)
    colnames(bw_df) <- c("name", "sample", "count")
    bw_df$window_type <- factor(rep("bait_windows", nrow(bw_df)))
    nbw_df <- reshape2::melt(nbw)
    colnames(nbw_df) <- c("name", "sample", "count")
    nbw_df$window_type <- factor(rep("non_bait_windows", nrow(nbw_df)))
    df <- rbind(bw_df,nbw_df)
    df <- tidyr::separate(df, name, c("chromosome", "start", "stop"), sep="[:-]" )
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

plot.atacr <- function(atcr){

#histogram of coverages by sample and window type
p1 <- atacr::coverage_histogram(atcr)

#density of coverage by chromosome region, bait windows
p2 <-  atacr::chromosome_coverage(atcr)


return(gridExtra::grid.arrange(p1,p2, nrow=2))



}
