#' @importFrom stats rlnorm

#' simulate counts and return an atacr object
#' @export
simulate_counts <- function() {
  # Each of the six sets of counts follows a mixed distribution of 10 counts drawn from a log-normal distribution with logmean 4 and SD 1, and 40 counts with logmean 10 and SD 1. This mimics the enrichment pattern we see with capture enriched data. 10 of the counts are multiplied by a value drawn from the normal distribution with mean 2 and SD 1 so can appear differentially expressed. These counts represent bait-windows - regions of the genome for which baits were designed and selected.
  num_windows = 100 #50 bait windows, 50 non bait windows
  reps = 3


  col_data <- S4Vectors::DataFrame(Treatment = c(rep("control", reps), rep("treatment", reps)),
    col.names = c(
      sprintf("control_%03d", 1:reps),
      sprintf("treatment_%03d", 1:reps)
    ))


  row_ranges <- GenomicRanges::GRanges(
    rep("synth_chrom", num_windows),
    IRanges::IRanges(seq(1, (num_windows * 50) , by = 50), width = 50),
    strand = sample(c("+", "-"), num_windows, TRUE),
    feature_id = sprintf("window_%06d", 1:num_windows)
  )

  names(row_ranges) <- sprintf("window_%06d", 1:num_windows)

  a <-
    floor(c(
      rlnorm(10, meanlog = 4, sdlog = 1),
      rlnorm(40, meanlog = 10, sdlog = 1)
    )) #basic two peak dist
  b <- floor(a * abs(rnorm(50, 1, sd = 1)))
  c <- floor(a * abs(rnorm(50, 1, sd = 1)))
  d <-
    floor(a * abs(c(rnorm(10, 2, sd = 1), rnorm(40, 1, sd = 1))))
  e <-
    floor(a * abs(c(rnorm(10, 2, sd = 1), rnorm(40, 1, sd = 1))))
  f <-
    floor(a * abs(c(rnorm(10, 2, sd = 1), rnorm(40, 1, sd = 1))))

  blank <- rep(0, 50)
  a <- c(a, blank)
  b <- c(b, blank)
  c <- c(c, blank)
  d <- c(d, blank)
  e <- c(e, blank)
  f <- c(f, blank)

  counts <-
    data.frame(
      control_001 = a,
      control_002 = b,
      control_003 = c,
      treatment_001 = d,
      treatment_002 = e,
      treatment_003 = f
    )

  counts <- as.matrix(counts[sample(nrow(counts)), ])
  row.names(counts) <- sprintf("window_%06d", 1:num_windows)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    rowRanges = row_ranges,
    colData = col_data
  )

  r <- list()
  class(r) <- c("atacr", "list")
  r$whole_genome <- se
  r$treatments <- c(rep("control", 3), rep("treatment", 3))
  r$sample_names <-
    c(sprintf("control_%03d", 1:3),
      sprintf("treatment_%03d", 1:3))
  r$bam_files <- "no.bam"

  bw <- which(counts[, 'control_001'] > 0)

  start_pos <- (bw - 1) * 50

  end_pos <- start_pos + 49
  seq_names <-  rep("synth_chrom", length(bw))
  r$bait_regions <- GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(seq_names),
    ranges = IRanges::IRanges(
      start_pos,
      end = end_pos,
      names = sprintf("bait_%02d", 1:length(bw))
    )
  )


  r$bait_windows <- r$whole_genome[bw,]
  r$non_bait_windows <- r$whole_genome[!bw,]
  r$whole_genome@rowRanges@ranges@NAMES <-
    as.character(r$whole_genome@rowRanges)
  r$bait_windows@rowRanges@ranges@NAMES <-
    as.character(r$bait_windows@rowRanges)
  r$non_bait_windows@rowRanges@ranges@NAMES <-
    as.character(r$non_bait_windows@rowRanges)
  r$dataframe <- as.data.frame(r)
  return(r)
}
