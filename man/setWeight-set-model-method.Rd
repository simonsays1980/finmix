% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model.R
\name{setWeight<-,model-method}
\alias{setWeight<-,model-method}
\title{Setter method of \code{model} class.}
\usage{
\S4method{setWeight}{model}(object) <- value
}
\arguments{
\item{object}{An \code{model} object.}

\item{value}{An matrix specifying the weights.}
}
\value{
The \code{model} object with slot \code{weight} set to \code{value}.
}
\description{
Sets a value for the \code{weight} slot.
}
\examples{
# Generate an default mixture model.
f_model <- model()
# Set the slot.
setWeight(f_model) <- matrix(c(0.4, 0.6), nrow = 1)

}
\seealso{
\itemize{
\item \linkS4class{model} for the class definition
}
}
\keyword{internal}
