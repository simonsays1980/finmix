% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prior.R
\name{getType,prior-method}
\alias{getType,prior-method}
\title{Getter method of \code{prior} class.}
\usage{
\S4method{getType}{prior}(object)
}
\arguments{
\item{object}{An \code{prior} object.}
}
\value{
The \code{type} slot of the \code{object}.
}
\description{
Returns the \code{type} slot.
}
\examples{
# Generate a prior object. 
f_prior <- prior()
# Get the slot.
getType(f_prior)
}
\keyword{internal}
