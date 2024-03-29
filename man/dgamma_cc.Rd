% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{dgamma_cc}
\alias{dgamma_cc}
\title{Computes the density of the Gamma distribution}
\usage{
dgamma_cc(values, shape, rate)
}
\arguments{
\item{values}{A matrix of dimension \verb{M x K} for which the density
should be calculated.}

\item{shape}{A vector of dimension \verb{K x 1} with Gamma shape parameters.}

\item{rate}{A vector of dimension \verb{K x 1} with Gamma rate parameters.}
}
\value{
A matrix of Gamma density values for each pair of parameters
in a column.
}
\description{
For each shape and rate parameter pair the gamma density is computed.
Inside the function the unsafe access functions of Armadillo \code{at()} and
\code{unsafe_col()} are used, so now boundary check is performed. In each step
the \code{lngamma()} function from Rcpp's \code{R} namespace is used. At this time
unused.
}
