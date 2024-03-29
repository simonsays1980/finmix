% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prior.R
\name{getPar,prior-method}
\alias{getPar,prior-method}
\title{Getter method of \code{prior} class.}
\usage{
\S4method{getPar}{prior}(object)
}
\arguments{
\item{object}{An \code{prior} object.}
}
\value{
The \code{par} slot of the \code{object}.
}
\description{
Returns the \code{par} slot.
}
\examples{
# Generate a prior object. 
f_prior <- prior()
# Get the slot.
getPar(f_prior)
}
\keyword{internal}
