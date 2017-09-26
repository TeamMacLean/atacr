#' estimates Goodness of Fit for each row in a count matrix
#' See:
#' https://haroldpimentel.wordpress.com/?s=TMM#paperList
#' https://academic.oup.com/biostatistics/article/13/3/523/248016/Normalization-testing-and-false-discovery-rate
#' https://github.com/cran/PoissonSeq/blob/3d9bc4b1744cb45714d4442b5a879b6e0c68b4a2/R/ps_other.R
#' @param mat a count matrix usually from SummarizedExperiment::assay()
#' @return a named vector of GoF estimates
gof <- function(mat){

  pseudo_val <- 1e-10

  shats <- colSums(mat) / sum(mat)
  x_shats <- rowSums(mat) %*% t(shats)
  gof <- rowSums((mat - x_shats) ^ 2 / (x_shats + pseudo_val))
  return(gof)
}

#' estimates Goodness of Fit from atacr object
#' @export
#' @param atacr a list of SummarizedExperiment objects from atacr::make_counts()
#' @param which the subdivision of the genome to calculate GoF either 'whole_genome', 'bait_windows' or 'non_bait_windows'
#' @return the original atacr object with a new slot - 'gofs' - a named vector of each windows GoF estimate.
estimate_GoFs <- function(atacr, which = "bait_windows"){
  mat <- SummarizedExperiment::assay(atacr[[which]])
  atacr$gofs <- gof(mat)
  return(atacr)
}

#' Depth estimation, directly from https://github.com/cran/PoissonSeq/blob/master/R/ps_cmeans.R
#' @param n a matrix
#' @param iter, runs of the Depth finder.
#' @return list of depths and means
Est.Depth <- function(n, iter=5)
{
  SMALL.VAL <- 1e-8
  cmeans <- colSums(n) / sum(n)
  keep <- NULL

  for (i in 1 : iter)
  {
    n0 <- rowSums(n) %*% t(cmeans)
    prop <- rowSums((n - n0) ^ 2 / (n0 + SMALL.VAL))
    qs <- quantile(prop, c(0.25, 0.75))
    keep <- (prop >= qs[1]) & (prop <= qs[2])

    cmeans <- colMeans(n[keep, ])
    cmeans <- cmeans / sum(cmeans)
  }

  return(list(cmeans=cmeans, keep=keep))
}



#' estimates sequencing depths based on windows with smallest GoF
#' @export
#' @param atacr a list of SummarizedExperiment objects from atacr::make_counts()
#' @param which the subdivision of the genome to calculate GoF either 'whole_genome', 'bait_windows' or 'non_bait_windows'
#' @return - a named vector of each windows GoF estimate.
get_GoF_factors <- function(atacr, which = "bait_windows"){

  m <- SummarizedExperiment::assay(atacr[[which]])
  seq.depth <- Est.Depth(n=m, iter=5)$cmeans
  seq.depth <- 1 / (exp(log(seq.depth) - mean(log(seq.depth))) )
  return(seq.depth)
}

#' find control windows by convergence method in https://academic.oup.com/biostatistics/article/13/3/523/248016/Normalization-testing-and-false-discovery-rate
#' @export
#' @param atacr a list of SummarizedExperiment objects from atacr::make_counts()
#' @param which the subdivision of the genome to calculate GoF either 'whole_genome', 'bait_windows' or 'non_bait_windows'
#' @return a character vector of window names
find_controls_by_GoF <- function(atacr, which = "bait_windows"){

  m <- SummarizedExperiment::assay(atacr[[which]])
  controls <- rownames(m[Est.Depth(n = m, iter = 5)$keep,])
  return(controls)

}


#' performs a whole library size normalisation of the selected set of windows, calculates a median virtual experiment and normalises to that
#' @export
#' @param data a list of SummarizedExperiment objects from atacr::make_counts()
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

average_matrix_by_sample <- function(data, which = "bait_windows") {
  l <- list()
  m <- SummarizedExperiment::assay(data[[which]])
  for (treatment in unique(data$treatments) ){
    samples <- names_from_treatment(data, treatment)
    l[[treatment]] <- apply(m[,samples], 1, mean)
  }
  return(do.call(cbind, l))

}


names_from_treatment <- function(data, treatment){
  return(data$sample_names[which(data$treatments == treatment)] )
}

#' return list of treatment names
#' @export
#' @param data an atacr object
#' @return char vector of unique treatment names
treatments <- function(data){
  return( unique(data$treatments))
}

treatment_from_name <- function(data, sample_name){
  return(data$treatments[which(data$sample_names == sample_name)] )
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
library_size_scaling_factors <- function( se ){ # nocov start
  sample_matrix <- SummarizedExperiment::assay( se )
  return(get_scaling_factors(sample_matrix))
} # nocov end

get_scaling_factors <- function( sample_matrix ){
  mve_sum <- sum(median_virtual_experiment( sample_matrix ))
  scaling_factors <- sapply(colSums( sample_matrix ), function(x){ mve_sum / x })
  return(scaling_factors)
}

scale_normalise <- function( sample_matrix, scaling_factors){ #nocov start
  scaled <- sample_matrix %*% diag(scaling_factors)
  return(scaled)
} #nocov end

#' normalise by a provided set of scaling factors
#' @export
#' @param data a list of SummarizedExperiment objects from atacr::make_counts()
#' @param which the subdivision of the genome to calculate correlations either 'whole_genome', 'bait_windows' or 'non_bait_windows'
#' @param scaling_factors a vector of scaling factors to normalise by
#' @return a SummarizedExperiment with scale normalised window values
scale_factor_normalise <- function(data, which = "bait_windows", scaling_factors = NULL){
  se <-  data[[which]]
  normalised_sample_matrix <- scale_normalise(SummarizedExperiment::assay( se ), scaling_factors)
  se_copy <- se
  SummarizedExperiment::assay(se_copy) <- normalised_sample_matrix
  return(se_copy)
}


#' extract scaling factors from control windows (often from a file of control gene positions)
#' @export
#' @param se a SummarizedExperiment object
#' @param window_file a text file containing the positions of control window/gene ranges
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
#' @param window_file a text file containing the positions of control window/gene ranges
#' @return SummarizedExperiment object, a copy of se with normalised values
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

#' normalise counts by window width (counts / window width)
#' @export
#' @param data a list of SummarizedExperiment objects from atacr::make_counts()
#' @param which the subset of the data to normalise. Default = bait_windows
#' @return SummarizedExperiment object with normalised counts
normalise_by_window_width <- function(data, which = "bait_windows"){
  widths <- data[[which]]@rowRanges@ranges@width
  return(
    SummarizedExperiment::assay(data[[which]]) / widths
  )
}
