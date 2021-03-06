% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/differentials.R
\name{estimate_bayes_factor}
\alias{estimate_bayes_factor}
\title{Estimate Bayes Factor and significantly different windows}
\usage{
estimate_bayes_factor(atacr, treatment_a, treatment_b, which = "bait_windows",
  factor = 4)
}
\arguments{
\item{atacr}{an atacr object}

\item{treatment_a}{the first treatment to consider}

\item{treatment_b}{the second treatment to consider}

\item{which}{the subset of windows to consider}

\item{factor}{the BayesFactor at which to mark window as significant}
}
\value{
a dataframe
}
\description{
Estimate Bayes Factor and significantly different windows
}
