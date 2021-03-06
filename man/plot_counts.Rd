% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/figures.R
\name{plot_counts}
\alias{plot_counts}
\title{Plot distribution of counts in given data set}
\usage{
plot_counts(data, which = "bait_windows", log10 = TRUE)
}
\arguments{
\item{data}{a list of SummarizedExperiment objects from atacr::make_counts()}

\item{which}{the subdivision of the genome to plot}

\item{log10}{log 10 the counts for plotting.}
}
\value{
ggplot2 plot
}
\description{
Plot distribution of counts in given data set
}
