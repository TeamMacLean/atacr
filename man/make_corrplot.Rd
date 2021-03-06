% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/figures.R
\name{make_corrplot}
\alias{make_corrplot}
\title{generate corrplot from matrix of counts}
\usage{
make_corrplot(counts, method = "pearson")
}
\arguments{
\item{counts}{a matrix of counts}

\item{method}{the correlation method to use, any supported by `cor()` is useable}
}
\value{
ggplot2 plot
}
\description{
generate corrplot from matrix of counts
}
