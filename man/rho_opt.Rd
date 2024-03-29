% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{rho_opt}
\alias{rho_opt}
\title{Quasi optimal \eqn{\rho} function}
\usage{
rho_opt(x, c)
}
\arguments{
\item{x}{numeric vector with positive values.}

\item{c}{tunning constant.}
}
\value{
Numeric vector with quasi optimal \eqn{\rho} computation for each element
of x.
}
\description{
Quasi optimal \eqn{\rho} function
}
\references{
[1] Salibian-Barrera, M., Willems, G., & Zamar, R. (2008). The fast-tau
estimator for regression. Journal of Computational and GraphicalStatistics,
17(3), 659-682.
}
