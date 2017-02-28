#' performs a whole library size normalisation of the selected set of windows, calculates a median virtual experiment and normalises to that
#'  @export
#'  @param data a list of SummarizedExperiment objects from atacr::make_counts()
#' @param which the subdivision of the genome to calculate correlations either 'whole_genome', 'bait_windows' or 'non_bait_windows'
#' @param by_treatment (FALSE) will group the assay into different treatments and normalise each separately - assumes that within treatment groups the samples should show little difference, but between sample treatment groups could show lots of difference between windows.
#' @return a SummarizedExperiment object with a new, normalised assay matrix
library_size_normalisation <- function(data, which = "bait_windows", by_treatment = FALSE ){
    l <- list()
    d <- data[[which]]
    if ( !by_treatment ){
      return( library_size_normalisation_internal( data[[which]] ) )
    }
    else {
      for (treatment in unique(data$treatments) ){
        samples <- names_from_treatment(data, treatment)
        treat_norm <- library_size_normalisation_internal( data[[which]][,samples] )
        l[[treatment]] <- SummarizedExperiment::assay( treat_norm )
      }
        full_mat <- do.call(cbind, l)
        SummarizedExperiment::assay(d) <- full_mat
        return( d )
    }


}

names_from_treatment <- function(data, treatment){
  return(data$sample_names[which(data$treatments == treatment)] )
}
#' do a library size normalisation
#' @param se a SummarizedExperiment object such as 'bait_windows' from atacr::make_counts()
library_size_normalisation_internal <- function(se){

  scaling_factors <- library_size_scaling_factors( se )
  normalised_sample_matrix <- scale_normalise(SummarizedExperiment::assay( se ), scaling_factors)
  se_copy <- se
  SummarizedExperiment::assay(se_copy) <- normalised_sample_matrix
  return(se_copy)
}

#' calculate scaling factors for library size
#' @export
#' @param se a SummarizedExperiment object such as 'bait_windows' from atacr::make_counts()
library_size_scaling_factors <- function( se ){
  sample_matrix <- SummarizedExperiment::assay( se )
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
#' @param by_treatment (FALSE) will group the assay into different treatments and normalise each separately - assumes that within treatment groups the samples should show little difference, but between sample treatment groups could show lots of difference between windows.
#' @return a vector of scaling factors from control genes
control_window_scaling_factors <- function( se, window_file){
  control_window_regions <- get_bait_regions_from_text( window_file )
  keep <- IRanges::overlapsAny( SummarizedExperiment::rowRanges( se ), control_window_regions )
  control_windows <- se[keep, ]
  sample_matrix <- SummarizedExperiment::assay( control_windows )
  return(get_scaling_factors(sample_matrix))

}
#' do a control window scaling normalisation
#' @param se a SummarizedExperiment object such as 'bait_windows' from atacr::make_counts()
control_window_normalise_internal <- function( se, window_file ){
  scaling_factors <- control_window_scaling_factors( se, window_file)
  normalised_sample_matrix <- scale_normalise(SummarizedExperiment::assay( se ), scaling_factors)
  se_copy <- se
  SummarizedExperiment::assay(se_copy) <- normalised_sample_matrix
  return(se_copy)

}
#' performs control window based normalisation of the selected set of windows, calculates a median virtual experiment and normalises to that
#' @export
#' @param data a list of SummarizedExperiment objects from atacr::make_counts()
#' @param which the subdivision of the genome to calculate correlations either 'whole_genome', 'bait_windows' or 'non_bait_windows'
#' @param window_file a text file containing the positions of control window/gene ranges
#' @param by_treatment should normalisation be done by all experiments (one median virtualexperiment to compare all samples to) OR should normalisation be done by each treatment type (one median virtual experiment for each different treatment type)
#' @return a vector of scaling factors from control genes
control_window_normalise <- function(data, window_file, which = "bait_windows", by_treatment = FALSE ){
    d <- data[[which]]
    l <- list()
    if(!by_treatment){
      return( control_window_normalise_internal(data[[which]], window_file ))
    }
    else{
      for (treatment in unique(data$treatments) ){
        samples <- names_from_treatment(data, treatment)
        treat_norm <- control_window_normalise_internal( data[[which]][,samples], window_file )
        l[[treatment]] <- SummarizedExperiment::assay( treat_norm )
      }
      full_mat <- do.call(cbind, l)
      SummarizedExperiment::assay(d) <- full_mat
      return( d )
    }
}
