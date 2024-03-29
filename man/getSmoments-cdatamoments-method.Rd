% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cdatamoments.R
\name{getSmoments,cdatamoments-method}
\alias{getSmoments,cdatamoments-method}
\title{Getter method of \code{cdatamoments} class.}
\usage{
\S4method{getSmoments}{cdatamoments}(object)
}
\arguments{
\item{object}{An \code{cdatamoments} object.}
}
\value{
The \code{smoments} slot of the \code{object}.
}
\description{
Returns the \code{smoments} slot.
}
\examples{
# Generate an exponential mixture model with two components.
f_model <- model("exponential", par = list(lambda = c(0.3, 0.7)), K = 2)
# Simulate data from the model.
f_data <- simulate(f_model)
# Calculate the mixture moments.
f_datamoms <- datamoments(f_data)
# Get the moments for the included indicators of the data. 
getSmoments(f_datamoms)

}
\seealso{
\itemize{
\item \linkS4class{datamoments} for the base class for model moments
\item \code{\link[=datamoments]{datamoments()}} for the constructor of the \code{datamoments} class family
}
}
\keyword{internal}
