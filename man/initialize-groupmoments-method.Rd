% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/groupmoments.R
\name{initialize,groupmoments-method}
\alias{initialize,groupmoments-method}
\title{Initializer of the \code{groupmoments} class}
\usage{
\S4method{initialize}{groupmoments}(.Object, ..., value)
}
\arguments{
\item{.Object}{An object: see the "initialize Methods" section in
\link{initialize}.}

\item{...}{Arguments to specify properties of the new object, to be passed
to \code{initialize()}.}

\item{model}{A finmix \linkS4class{fdata} object containing the observations.}
}
\description{
Only used implicitly. The initializer calls a function \code{generateMoments()}
object. to generate in the initialization step the moments for a passed-in
\code{fdata} object.
}
\seealso{
\itemize{
\item \link{Classes_Details} for details of class definitions, and
\item \link{setOldClass} for the relation to S3 classes
}
}
\keyword{internal}
