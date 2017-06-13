#' generate a synthetic count matrix
#' @param num_windows the number of windows to generate (rows in the matrix)
#' @param reps the number of control and treatmentsexperiments to generate
#' @param diff_by the fold-change of the differentially accessible windows
#' @param prop_baits the proportion of the windows that are used as baits
#' @param diff_prop proportion of windows to make differentially accesible
#' @param expected_coverage expected count per window (mu of neg binomial)
#' @param dispersion_factor expected increase in size of dispersion relative to mean
#' @return a RangedSummarizedExperiment object with count matrix as specified
#' @export
make_synthetic_count_SE <- function(num_windows, prop_baits = 0.5, diff_prop=0.1, reps = 3,  diff_by = 2, expected_coverage = 30, dispersion_factor = 100){ # nocov start

  #the number of rows in the bait matrix
  num_bait_rows <- ceiling(num_windows * prop_baits)

  #the number of data points in the bait matrix
  num_bait_data_points <- num_bait_rows * reps

  #the number of data points in the non_bait matrix
  num_non_bait_data_points = (num_windows * reps) - num_bait_data_points

  #the number of rows in the non bait matrix
  num_non_bait_rows <- (num_windows * (1 - prop_baits))

    #set up probs of 0,1,2,3 and 4 counts in non-bait windows
  non_bait_probs = c(rep(0, 13), rep(1,3), rep(2,2), rep(3,1), rep(4,1))

  #generates non_bait window counts
  j <- matrix(sample(non_bait_probs,
                     num_non_bait_data_points,
                     replace = TRUE),
              nrow = num_non_bait_rows)
  k <- matrix(sample(non_bait_probs,
                     num_non_bait_data_points,
                     replace = TRUE),
              nrow = num_non_bait_rows)
  non_bait_counts <- cbind(j,k)

  m <- matrix(round(rnbinom( num_bait_data_points,
                       mu=expected_coverage,
                       size = (expected_coverage * dispersion_factor)
                     )),
               nrow = num_bait_rows  )

  n <- m

  #multiply everything by a random co-efficient to simulate experimental noise
  n <- round(n * runif(n, 0.8,1.2))

  # set the first diff_prop windows to a random number with mean (current value * diff_by)
  if (diff_prop != 0){
    num_da <- ceiling( (num_windows * prop_baits) * diff_prop )
    n[1:num_da,] <- round(n[1:num_da,] * rnorm(num_da, mean=diff_by, sd=diff_by/2 ))
  }

  bait_counts <- cbind(m,n)
  bait_counts[ bait_counts < 0 ] <- 0

  counts <- rbind(bait_counts, non_bait_counts)

  col_data <- S4Vectors::DataFrame(
    Treatment=c(rep("control", reps),rep("treatment", reps) ),
    row.names = c( sprintf("control_%03d", 1:reps), sprintf("treatment_%03d", 1:reps) )
  )

  row_ranges <- GenomicRanges::GRanges(
    rep("synth_chrom", num_windows),
    IRanges::IRanges(seq(1, (num_windows  * 50) , by=50), width=50),
    strand = sample(c("+", "-"), num_windows, TRUE),
    feature_id = sprintf("window_%06d", 1:num_windows)
  )
  names(row_ranges) <- sprintf("window_%06d", 1:num_windows)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts=counts),
    rowRanges = row_ranges,
    colData = col_data
  )
  return(se)
}

label_bait_regions <- function(num_windows, prop_baits){

  num_baits <-  ceiling(num_windows * prop_baits)
  bait_names <-  sprintf("bait_%02d", 1:num_baits)
  start_pos <-  seq(1,(num_baits * 50), by=50)
  end_pos <-  start_pos + 49
  seq_names <-  rep("synth_chrom",num_baits)
  bait_regions <- GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(seq_names),
    ranges = IRanges::IRanges(start_pos, end = end_pos, names = bait_names)
  )
  x <- 1
  return(bait_regions)
}

make_synthetic_data <- function(num_windows, prop_baits = 0.5, diff_prop=0.1, reps = 3,  diff_by = 2, expected_coverage = 30, dispersion_factor = 100){
  d <- atacr::make_synthetic_count_SE(num_windows)
  r <- list()
  class(r) <- c("atacr", "list")
  r$whole_genome <- d
  r$treatments <- c(rep("control", 3), rep("treatment", 3))
  r$sample_names <- c( sprintf("control_%03d", 1:3), sprintf("treatment_%03d", 1:3) )
  r$bam_files <- "no.bam"


  ## name samples
  colnames(r$whole_genome) <- r$sample_names
  r$bait_regions <- label_bait_regions(num_windows, prop_baits)


  keep <- IRanges::overlapsAny( SummarizedExperiment::rowRanges( r$whole_genome ), r$bait_regions )
  r$bait_windows <- r$whole_genome[keep, ]
  r$non_bait_windows <- r$whole_genome[!keep, ]
  r$whole_genome@rowRanges@ranges@NAMES <- as.character(r$whole_genome@rowRanges)
  r$bait_windows@rowRanges@ranges@NAMES <- as.character(r$bait_windows@rowRanges)
  r$non_bait_windows@rowRanges@ranges@NAMES <- as.character(r$non_bait_windows@rowRanges)
  r$dataframe <- as.data.frame(r)
  return(r)
} #nocov end

