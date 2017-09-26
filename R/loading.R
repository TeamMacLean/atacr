## loads in the bams and the list of bait windows.
## uses the read size to make valid windows
## works out the windows and percent of reads outside the windows
## works out coverage stats per experiment

#' load BAM files and calculate window coverage
#' @export
#' @param window_file A filename of a CSV file with the bait regions
#' @param sample_treatment_file A filename of a CSV file that lists treatments, samples and bam file paths
#' @param width an integer of the width of the bins the bait regions will be divided into
#' @param filter_params a params object from atacr::make_params()  that define how reads will be extracted from the BAM files. Optionally, for greater control, either a csaw::readParam() (for ATACseq) or Rsamtools::ScanBamParam() object for RNASeq can be provided. See http://bioconductor.org/packages/release/bioc/manuals/csaw/man/csaw.pdf or https://www.rdocumentation.org/packages/Rsamtools/versions/1.24.0/topics/ScanBamParam for details
#' @param is_rnaseq a boolean stating whether this is RNASeq data. Default = FALSE
#' @param gene_id_col a character string stating which attribute name to take from the final column of the GFF file to use for the window name in RNASeq data. Usually this is the name of the gene. Default = ID.
#' @param with_df attach a dataframe version of the data Default = FALSE
#' @return a list of metadata and RangedSummarizedExperiment objects with read count in windows for whole genome, bait windows and non-bait windows for each sample
make_counts <-
  function(window_file,
    sample_treatment_file,
    width = 50,
    filter_params = make_params(), #csaw::readParam(minq = 50),
    with_df = FALSE,
    is_rnaseq = FALSE,
    gene_id_col = "ID") {
    result <- list()
    class(result) <- c("atacr", "list")

    sample_treatment_file_mapping <-
      read_experiment_info(sample_treatment_file)
    result$bam_files <-
      as.character(sample_treatment_file_mapping$bam_file_path)
    result$treatments <-
      as.character(sample_treatment_file_mapping$treatment)
    result$sample_names <-
      as.character(sample_treatment_file_mapping$sample_name)

    if (!is_rnaseq) {
      result <- load_atac(result, width, filter_params, window_file)
    }
    else {
      result <- load_rnaseq(result, filter_params, window_file)
    }

    if (with_df) {
      result$dataframe <- as.data.frame(result)
    }
    return(result)
  }

#' set read filters for counting from the BAM file.
#' @export
#' @param paired_map Should reads only be included if they are aligned in pairs. Default = TRUE
#' @param minq The minimum mapping quality to retain a read. Default = 20
#' @param dedup Should removal of PCR duplicates be performed. Default = TRUE
#' @return a named vector of class "atacr_params"
make_params <- function(paired_map = TRUE, minq = 30, dedup = TRUE){

  params <- c(paired_map, minq, dedup)
  names(params) <- c("paired_map", "minq", "dedup")
  class(params) <- c("atacr_params")
  return(params)

}

#' reads a csv file containing the bait regions
#' @param file_name path to a csv file containing the bait regions. File must have a header with columns `bait_name`, `seq_name`, `start`, `end`.
#' @return GenomicRanges object of bait regions
get_bait_regions_from_text <- function(file_name) {
  df <- read.csv(file_name, sep = ",", header = TRUE)

  bait_regions <- GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(df$seq_name),
    ranges = IRanges::IRanges(
      df$start_pos,
      end = df$end_pos,
      names = df$bait_name
    )
  )

  return(bait_regions)

}

#' reads a gff file containing the bait regions
#' @param file_name path to the file containing the bait regions
#' @return GenomicRanges object of bait regions
get_bait_regions_from_gff <- function(file_name) {
  gff <- rtracklayer::import.gff(file_name)
  bait_regions <- as(gff, "GRanges")
  bait_regions <- bait_regions[bait_regions$type %in% c("gene")]
  return(bait_regions)
}

#' format a csaw::readParam object from the atacr::make_params() object
#' @param p an object returned from atacr::make_params()
#' @return a csaw::readParam object
make_csaw_params <- function(p){
  return(
    csaw::readParam(
      minq = p["minq"],
      dedup = p["dedup"],
      pe = ifelse(p["paired_map"], "both", "none")
    )
  )
}

#' populate the result object with the RangedSummarizedExperiment from the bam files from ATAC seq data. Called from make_counts() when is_rnaseq == FALSE.
#' @param result list from make_counts()
#' @param width an integer of the width of the bins the bait regions will be divided into
#' @param filter_params a params object, described in atacr::make_counts()
#' @param window_file  a filename of a CSV file with the bait regions
#' @return a list with window counts for bait/non-bait windows
load_atac <- function(result, width, filter_params, window_file) {

  if ("atacr_params" %in% class(filter_params) ) {
    filter_params <- make_csaw_params(filter_params)
  }

  result$whole_genome <-
    csaw::windowCounts(
      result$bam_files,
      bin = TRUE,
      filter = 0,
      width = width,
      param = filter_params
    )

  ## name samples
  colnames(result$whole_genome) <- result$sample_names

  ## collect bait and non bait regions
  result$bait_regions <- get_bait_regions_from_text(window_file)

  keep <-
    IRanges::overlapsAny(SummarizedExperiment::rowRanges(result$whole_genome),
      result$bait_regions)

  result$bait_windows <- result$whole_genome[keep,]
  result$non_bait_windows <- result$whole_genome[!keep,]

  ## name the windows
  result$whole_genome@rowRanges@ranges@NAMES <-
    as.character(result$whole_genome@rowRanges)
  result$bait_windows@rowRanges@ranges@NAMES <-
    as.character(result$bait_windows@rowRanges)
  result$non_bait_windows@rowRanges@ranges@NAMES <-
    as.character(result$non_bait_windows@rowRanges)

  return(result)
}

#' format a rsamtools::scanBam object from the atacr::make_params() object
#' @param p an object returned from atacr::make_params()
#' @param example_bam a filename pointing to a BAM file from which genome size can be taken
#' @return an rsamtools::scanBamParam object
make_scanBamParam <- function(p, example_bam){

  seqnames <- seqlength <- NULL

  ranges <- Rsamtools::idxstatsBam(example_bam)
  ranges <- dplyr::mutate(ranges, start = 1)
  ranges <- dplyr::rename(ranges, seqname = seqnames, end = seqlength)
  ranges <- unlist(GenomicRanges::makeGRangesListFromDataFrame(ranges))

  return(Rsamtools::ScanBamParam(
    flag = Rsamtools::scanBamFlag(
      isDuplicate = !p["dedup"],
      isProperPair = p["paired_map"]
    ),
    mapqFilter = p["minq"],
    which = ranges
  ))
}

#' populate the result object with the RangedSummarizedExperiment from the bam files from RNA seq data. Called from make_counts() when is_rnaseq == TRUE.
#' @param result list from make_counts()
#' @param filter_params a params object, described in atacr::make_counts()
#' @param window_file  a filename of a CSV file with the bait regions
#' @param gene_id_col a character string stating which attribute name to take from the final column of the GFF file to use for the window name in RNASeq data. Usually this is the name of the gene. Default = ID.
load_rnaseq <-
  function(result,
    filter_params,
    window_file,
    gene_id_col = "ID") {

    if ("atacr_params" %in% class(filter_params) ) {
      filter_params <- make_scanBamParam(filter_params, result$bam_files[1])
    }

    bams <- Rsamtools::BamFileList(result$bam_files)
    names(bams) <- result$sample_names

    result$bait_regions <- get_bait_regions_from_gff(window_file)
    non_bait_regions <-
      GenomicRanges::gaps(result$bait_regions) #the intergene regions

    result$bait_windows <-
      GenomicAlignments::summarizeOverlaps(
        features = result$bait_regions,
        reads = bams,
        ignore.strand = T,
        param = filter_params
      )
    result$non_bait_windows <-
      GenomicAlignments::summarizeOverlaps(
        features = non_bait_regions,
        reads = bams,
        ignore.strand = T,
        param = filter_params
      )

    if (c(gene_id_col) %in% names(result$bait_regions@elementMetadata@listData)) {
      result$bait_windows@rowRanges@ranges@NAMES <-
        as.character(result$bait_regions@elementMetadata@listData[[gene_id_col]])
    }
    else {
      result$bait_windows@rowRanges@ranges@NAMES <- make_range_names(
        result$bait_regions@seqnames@values,
        result$bait_regions@ranges@start,
        result$bait_regions@ranges@width
      )
    }

    result$non_bait_windows@rowRanges@ranges@NAMES <- make_range_names(
      non_bait_regions@seqnames@values,
      non_bait_regions@ranges@start,
      non_bait_regions@ranges@width
    )

    result$bait_windows@rowRanges@elementMetadata@listData <- list()
    result$whole_genome <-
    rbind(result$bait_windows, result$non_bait_windows)
    colnames(result$whole_genome) <-
    colnames(result$bait_windows) <-
    colnames(result$non_bait_windows) <- result$sample_names

    return(result)
  }

make_range_names <- function(chr, start, width) {
  end <- start + width
  return(paste0(chr, ":", start, "-", end))
}

#' Loads in a CSV file describing treatment, samples and bam files
#' @param filename path and name of the file to load
read_experiment_info <- function(filename) {
  return(read.csv(filename, header = TRUE, sep = ","))
  ## read in file with columns 'treatment, sample_name, bam_file_path'
}

#' pulls lines out of a gff file based on identifierss provided
#' @param ids character vector of ids/names of feature to extract
#' @param gff path to gff file
#' @param type feature type of features to extract.
#' @param col column name of GFF file containing id to use (ID)
#' @param out_file path of file name to write. If NULL, no file is written
#' @param version which gff version to export (Default is "3")
#' @return GenomicRanges or NULL with GFF outfile.
extract_features_from_gff <- function(ids, gff, type = c("gene"), col="ID", out_file = NULL, version = "3"){
  gff <- rtracklayer::import.gff(gff, "GFF")
  gene_features <- as(gff, "GRanges")
  gene_features <- gene_features[ gene_features$type %in% type ]
  gene_features <- gene_features[ gene_features@elementMetadata@listData[[col]] %in% ids ]
  if ( is.null(out_file) ) {
    return(gene_features)
  }
  else{
    rtracklayer::export.gff(gene_features, out_file, version = version)
  }

}
