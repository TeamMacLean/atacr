% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/figures.R
\name{coverage_summary}
\alias{coverage_summary}
\title{Plot histograms of read counts by sample and window type}
\usage{
coverage_summary(data, which = NULL, sample = NULL, log_axis = TRUE)
}
\arguments{
\item{data}{a list of SummarizedExperiment objects from atacr::make_counts()}

\item{which}{the subdivision of the genome to plot (default = bait and non_bait windows)}

\item{sample}{the sample to plot (default = all )}

\item{log_axis}{use a log scale for the x-axis}
}
\value{
a ggplot2 object
}
\description{
Plot histograms of read counts by sample and window type
}
