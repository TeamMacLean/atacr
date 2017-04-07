## loads in the bams and the list of bait windows.
## uses the read size to make valid windows
## works out the windows and percent of reads outside the windows
## works out coverage stats per experiment

#' load BAM files and calculate window coverage
#' @export
#' @param window_file A filename of a CSV file with the bait regions
#' @param sample_treatment_file A filename of a CSV file that lists treatments, samples and bam file paths
#' @param width an integer of the width of the bins the bait regions will be divided into
#' @param filter_params a csaw::readParam() object that defines how reads will be extracted from the BAM files. See http://bioconductor.org/packages/release/bioc/manuals/csaw/man/csaw.pdf for details
#' @param with_df attach a dataframe version of the data Default = FALSE
#' @return a list of metadata and RangedSummarizedExperiment objects with read count in windows for whole genome, bait windows and non-bait windows for each sample
make_counts <- function(window_file, sample_treatment_file, width=50, filter_params = csaw::readParam(minq = 50), with_df = FALSE ){

  result <- list()
  class(result) <- c("atacr", "list")

  sample_treatment_file_mapping <-  read_experiment_info(sample_treatment_file)
  result$bam_files <- as.character(sample_treatment_file_mapping$bam_file_path)
  result$treatments <- as.character(sample_treatment_file_mapping$treatment)
  result$sample_names <- as.character(sample_treatment_file_mapping$sample_name)
  result$whole_genome <- csaw::windowCounts(result$bam_files, bin = TRUE, filter = 0, width = width, param = filter_params)


  ## name samples
  colnames(result$whole_genome) <- result$sample_names

  ## collect bait and non bait regions
  result$bait_regions <- get_bait_regions_from_text(window_file)

  keep <- IRanges::overlapsAny( SummarizedExperiment::rowRanges( result$whole_genome ), result$bait_regions )

  result$bait_windows <- result$whole_genome[keep, ]
  result$non_bait_windows <- result$whole_genome[!keep, ]

  ## name the windows
  result$whole_genome@rowRanges@ranges@NAMES <- as.character(result$whole_genome@rowRanges)
  result$bait_windows@rowRanges@ranges@NAMES <- as.character(result$bait_windows@rowRanges)
  result$non_bait_windows@rowRanges@ranges@NAMES <- as.character(result$non_bait_windows@rowRanges)

  if (with_df){
    result$dataframe <- as.data.frame(result)
  }
  return(result)
}

#' reads a text file containing the bait regions
get_bait_regions_from_text <- function(file_name){
  df <- read.csv(file_name, sep="\t", header=FALSE)
  df <- plyr::rename(df, c("V1"="bait_name", "V2"="seq_name", "V3"="col1", "V4"="col2"))
  df$start_pos <- apply(df[,3:4],1,min)
  df$end_pos <-  apply(df[,3:4],1,max)

  bait_regions <- GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(df$seq_name),
    ranges = IRanges::IRanges(df$start_pos, end = df$end_pos, names = df$bait_name)
  )

  return(bait_regions)

}




#' Loads in a CSV file describing treatment, samples and bam files
#' @param filename path and name of the file to load
read_experiment_info <- function(filename){
  return(read.csv(filename, header=TRUE, sep=",") )
  ## read in file with columns 'treatment, sample_name, bam_file_path'
}

