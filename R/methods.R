meta_summary <- function(atcr){
  samples = paste(unique(atcr$sample_names), collapse=",")
  treatments = paste(unique(atcr$treatments), collapse=",")
  sample_count = length(unique(atcr$sample_names))
  treat_count  = length(unique(atcr$treatments))
  return(cat(
    "ATAC-seq experiment of", treat_count, "treatments in", sample_count, "samples\n",
    "Treatments:", treatments, "\n",
    "Samples:", samples, "\n",
    "Windows used:", length(atcr$bait_regions), "\n"
    ))
}

print.atacr <- function(atcr){
  meta_summary(atcr)
}

summary.atacr <- function(atcr){
  meta <- meta_summary(atcr)
  on_target <- paste(capture.output(target_count_summary(atcr)), collapse="\n")
  coverage <- paste(capture.output(coverage_count_summary(atcr)), collapse="\n")
  return(cat(
    meta, "\n",
    "On/Off target read counts:\n",
    on_target,"\n",
    "Read depths:\n",
    coverage
    ) )


}

