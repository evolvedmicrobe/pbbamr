% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{toPhred}
\alias{toPhred}
\title{Function converts to Phred Scale, capping the error at a specified minimum}
\usage{
toPhred(errRates, minError = -Inf)
}
\arguments{
\item{minError}{The minimum allowed error rate.}

\item{errRate}{The error rate observed}
}
\value{
A vector of phred scaled error rates
}
\description{
Function converts to Phred Scale, capping the error at a specified minimum
}
\examples{
toPhred(10 ^ (-(1:4)), 1e-3)
}

