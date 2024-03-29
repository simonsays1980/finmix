% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model.R
\name{simulate,model-method}
\alias{simulate,model-method}
\title{Simulates data from a model.}
\usage{
\S4method{simulate}{model}(model, N = 100, varargin, seed = 0)
}
\arguments{
\item{model}{An S4 model object with specified parameters and weights.}

\item{N}{An integer specifying the number of values to be simulated.}

\item{varargin}{An S4 fdata object with specified variable dimensions, \code{r}
and repetitions \code{T}.}

\item{seed}{An integer specifying the seed for the RNG.}
}
\value{
An S4 fdata object holding the simulated values.
}
\description{
\code{simulate()} simulates values for a specified mixture model in an
S4 \code{model} object.
}
\examples{
\dontrun{
f_data <- simulate(model, 100)
}

}
\seealso{
\itemize{
\item \linkS4class{model} for the class definition
\item \linkS4class{fdata} for the class defining \code{finmix} data objects
}
}
\keyword{internal}
