% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{derpsi_opt}
\alias{derpsi_opt}
\title{derpsiOpt
the derivative of the psi function}
\usage{
derpsi_opt(t, c)
}
\arguments{
\item{t}{a numeric vector}

\item{c}{a tunning constant.}
}
\value{
rho'(\code{x}/ \code{c})
}
\description{
derpsiOpt
the derivative of the psi function
}
\examples{
val <- derpsi_opt(t = 0.5, c = 1.0)

}
