% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{wni}
\alias{wni}
\title{Weight function}
\usage{
wni(distances, c1, c2, s)
}
\arguments{
\item{distances}{numeric vector with distances from each point to its
cluster center.}

\item{c1}{tuning constant used for m scale estimation.}

\item{c2}{tuning constant used for \eqn{\tau} scale estimation.}

\item{s}{M scale for the same vector of distances.}
}
\value{
Numeric vector with the weight for each observation
}
\description{
Weight function
}
