#' performs a whole library size normalisation of the selected set of windows, calculates a median virtual experiment and normalises to that
#'  @export
#'  @param data a list of SummarizedExperiment objects from atacr::make_counts()
#' @param which the subdivision of the genome to calculate correlations either 'whole_genome', 'bait_windows' or 'non_bait_windows'
#' @return a SummarizedExperiment object with a new, normalised assay matrix
library_size_normalisation <- function(data, which = "bait_windows"){

  scaling_factors <- library_size_scaling_factors(data, which)
  normalised_sample_matrix <- scale_normalise(SummarizedExperiment::assay( data[[which]]), scaling_factors)
  se_copy <- data[[which]]
  SummarizedExperiment::assay(se_copy) <- normalised_sample_matrix
  return(se_copy)
}

#' calculate scaling factors for library size
#' @export
#' @param data a list of SummarizedExperiment objects from atacr::make_counts()
#' @param which the subdivision of the genome to calculate correlations either 'whole_genome', 'bait_windows' or 'non_bait_windows'
library_size_scaling_factors <- function(data, which = "bait_windows"){
  sample_matrix <- SummarizedExperiment::assay( data[[which]] )
  return(get_scaling_factors(sample_matrix))
}

get_scaling_factors <- function( sample_matrix ){
  mve_sum <- sum(atacr::median_virtual_experiment( sample_matrix ))
  scaling_factors <- sapply(colSums( sample_matrix ), function(x){ mve_sum / x })
  return(scaling_factors)
}

scale_normalise <- function( sample_matrix, scaling_factors){
  scaled <- sample_matrix %*% diag(scaling_factors)
  return(scaled)
}

#' extract scaling factors from control windows (often from a file of control gene positions)
#' @export
#' @param data a list of SummarizedExperiment objects from atacr::make_counts()
#' @param which the subdivision of the genome to calculate correlations either 'whole_genome', 'bait_windows' or 'non_bait_windows'
#' @param window_file a text file containing the positions of control window/gene ranges
#' @return a vector of scaling factors from control genes
control_window_scaling_factors <- function( data, window_file, which = "bait_windows" ){
  control_window_regions <- get_bait_regions_from_text( window_file )
  keep <- IRanges::overlapsAny( SummarizedExperiment::rowRanges( data[[which]] ), control_window_regions )

  control_windows <- data[[which]][keep, ]
  sample_matrix <- SummarizedExperiment::assay( control_windows )
  return(get_scaling_factors(sample_matrix))

}

#' performs control window based normalisation of the selected set of windows, calculates a median virtual experiment and normalises to that
#' @export
#' @param data a list of SummarizedExperiment objects from atacr::make_counts()
#' @param which the subdivision of the genome to calculate correlations either 'whole_genome', 'bait_windows' or 'non_bait_windows'
#' @param window_file a text file containing the positions of control window/gene ranges
#' @return a vector of scaling factors from control genes
control_window_normalise <- function(data, window_file, which = "bait_windows"){
  scaling_factors <- control_window_scaling_factors(data, window_file, which = which)
  normalised_sample_matrix <- scale_normalise(SummarizedExperiment::assay( data[[which]]), scaling_factors)
  se_copy <- data[[which]]
  SummarizedExperiment::assay(se_copy) <- normalised_sample_matrix
  return(se_copy)
}
