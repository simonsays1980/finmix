% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/binomialmodelmoments.R
\name{getExtrabinvar,binomialmodelmoments-method}
\alias{getExtrabinvar,binomialmodelmoments-method}
\title{Getter method of \code{binomialmodelmoments} class.}
\usage{
\S4method{getExtrabinvar}{binomialmodelmoments}(object)
}
\arguments{
\item{object}{An \code{binomialmodelmoments} object.}
}
\value{
The \code{extrabinvar} slot of the \code{object}.
}
\description{
Returns the \code{extrabinvar} slot.
}
\examples{
f_model <- model("binomial", par=list(p=c(0.3, 0.5)), 
                 weight=matrix(c(0.3, 0.7), nrow=1))
f_moments <- modelmoments(f_model)
getExtrabinvar(f_moments)

}
\seealso{
\itemize{
\item \linkS4class{modelmoments} for the base class for model moments
\item \code{\link[=modelmoments]{modelmoments()}} for the constructor of the \code{modelmoments} class family
}
}
\keyword{internal}
