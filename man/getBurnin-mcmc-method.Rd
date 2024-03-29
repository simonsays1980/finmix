% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc.R
\name{getBurnin,mcmc-method}
\alias{getBurnin,mcmc-method}
\title{Getter method of \code{mcmc} class.}
\usage{
\S4method{getBurnin}{mcmc}(object)
}
\arguments{
\item{object}{An \code{mcmc} object.}
}
\value{
The \code{burnin} slot of the \code{object}.
}
\description{
Returns the \code{burnin} slot.
}
\examples{
# Generate an mcmc object
f_mcmc <- mcmc()
# Get the slot
getBurnin(f_mcmc)

}
\seealso{
\itemize{
\item \linkS4class{mcmc} for the class definition
\item \code{\link[=mcmc]{mcmc()}} for the constructor of the \code{mcmc} class
}
}
\keyword{internal}
