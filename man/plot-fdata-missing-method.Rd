% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fdata.R
\name{plot,fdata,missing-method}
\alias{plot,fdata,missing-method}
\title{Plots the data}
\usage{
\S4method{plot}{fdata,missing}(x, y, dev = TRUE, ...)
}
\arguments{
\item{x}{An \code{fdata} object. Cannot be empty.}

\item{y}{Unused.}

\item{dev}{A logical indicating if the plot should be output via a graphical
device.}

\item{...}{Further arguments passed to the plotting functions \code{hist} or
\code{barplot}.}
}
\description{
\code{\link[=plot]{plot()}} plots the data in an \link{fdata} object by either a histogram in case of
continuous data or a barplot in case of discrete data.
}
\examples{
# Generate Poisson data and plot it. 
f_data <- fdata(y = rpois(100, 312), sim = TRUE)
plot(f_data)

}
\seealso{
\link{fdata} class
}
\keyword{internal}
