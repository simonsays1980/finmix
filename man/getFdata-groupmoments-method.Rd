% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/groupmoments.R
\name{getFdata,groupmoments-method}
\alias{getFdata,groupmoments-method}
\title{Getter method of \code{groupmoments} class.}
\usage{
\S4method{getFdata}{groupmoments}(object)
}
\arguments{
\item{object}{An \code{groupmoments} object.}
}
\value{
The \code{fdata} slot of the \code{object}.
}
\description{
Returns the \code{fdata} slot.
}
\examples{
# Generate a Poisson mixture model with two components.
f_model <- model("poisson", par = list(lambda = c(0.3, 0.7)), K = 2)
# Simulate data from the model.
f_data <- simulate(f_model)
# Calculate the mixture moments.
f_gmoments <- groupmoments(f_data)
# Get the data<
getFdata(f_gmoments)

}
\seealso{
\itemize{
\item \linkS4class{groupmoments} for the definition of the \code{groupmoments}
class
\item \code{\link[=groupmoments]{groupmoments()}} for the class constructor
}
}
\keyword{internal}
