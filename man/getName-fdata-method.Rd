% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fdata.R
\name{getName,fdata-method}
\alias{getName,fdata-method}
\title{Getter method of \code{fdata} class.}
\usage{
\S4method{getName}{fdata}(object)
}
\arguments{
\item{object}{An \code{fdata} object.}
}
\value{
The \code{name} slot of the \code{object}.
}
\description{
Returns the \code{name} slot of an \code{fdata} object.
}
\examples{
# Create an fdata object with Poisson data
f_data <- fdata(y = rpois(100, 312), sim = TRUE)
getName(f_data)

}
\seealso{
\link{fdata} for all slots of the \code{fdata} class
}
\keyword{internal}
