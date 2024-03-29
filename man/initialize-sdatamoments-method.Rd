% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sdatamoments.R
\name{initialize,sdatamoments-method}
\alias{initialize,sdatamoments-method}
\title{Initializer of the \code{sdatamoments} class}
\usage{
\S4method{initialize}{sdatamoments}(.Object, ..., value = fdata())
}
\arguments{
\item{.Object}{An object: see the "initialize Methods" section in
\link{initialize}.}

\item{...}{Arguments to specify properties of the new object, to be passed
to \code{initialize()}.}

\item{model}{A finmix \linkS4class{fdata} object containing the observations.}
}
\description{
Only used implicitly. The initializer calls the constructor for a
\linkS4class{groupmoments} object. to generate in the initialization
step the moments for a passed-in \code{fdata} object.
}
\seealso{
\itemize{
\item \link{Classes_Details} for details of class definitions, and
\item \link{setOldClass} for the relation to S3 classes
}
}
\keyword{internal}
