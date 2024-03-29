% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cdatamoments.R
\name{getHigher,cdatamoments-method}
\alias{getHigher,cdatamoments-method}
\title{Getter method of \code{cdatamoments} class.}
\usage{
\S4method{getHigher}{cdatamoments}(object)
}
\arguments{
\item{object}{An \code{cdatamoments} object.}
}
\value{
The \code{higher} slot of the \code{object}.
}
\description{
Returns the \code{higher} slot.
}
\examples{
# Generate an exponential mixture model with two components.
f_model <- model("exponential", par = list(lambda = c(0.3, 0.7)), K = 2)
# Simulate data from the model.
f_data <- simulate(f_model)
# Calculate the mixture moments.
f_datamoms <- datamoments(f_data)
# Use the getter.
getHigher(f_datamoms)

}
\seealso{
\itemize{
\item \linkS4class{datamoments} for the base class for model moments
\item \code{\link[=datamoments]{datamoments()}} for the constructor of the \code{datamoments} class family
}
}
\keyword{internal}
