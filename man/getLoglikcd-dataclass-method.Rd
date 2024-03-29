% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataclass.R
\name{getLoglikcd,dataclass-method}
\alias{getLoglikcd,dataclass-method}
\title{Getter method of \code{dataclass} class.}
\usage{
\S4method{getLoglikcd}{dataclass}(object)
}
\arguments{
\item{object}{An \code{dataclass} object.}
}
\value{
The \code{loglikcd} slot of the \code{object}.
}
\description{
Returns the \code{loglikcd} slot. Note that this slot is only non-null, if the
indicators have not been simulated.
}
\examples{
# Generate a Poisson mixture model with two components.
f_model <- model("poisson", par = list(lambda = c(0.3, 0.7)), K = 2)
# Simulate data from the model.
f_data <- simulate(f_model)
# Classify observations
f_dataclass <- dataclass(f_data, f_model, simS = FALSE)
getLoglikcd(f_dataclass)

}
\seealso{
\itemize{
\item \linkS4class{dataclass} for the base class
\item \code{\link[=dataclass]{dataclass()}} for the constructor of the \code{dataclass} class
}
}
\keyword{internal}
