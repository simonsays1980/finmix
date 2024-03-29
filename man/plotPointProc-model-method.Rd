% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model.R
\name{plotPointProc,model-method}
\alias{plotPointProc,model-method}
\title{Plots the point process of a finite model}
\usage{
\S4method{plotPointProc}{model}(x, dev = TRUE, ...)
}
\arguments{
\item{x}{An S4 model object with defined parameters and weights.}

\item{dev}{A logical indicating, if the plot should be shown in a graphical
device. Set to \code{FALSE}, if plotted to a file.}

\item{...}{Arguments to be passed to methods, such as graphical parameters
(see \link{par}).}
}
\value{
A scatter plot of weighted parameters.
}
\description{
\code{plotPointProc} plots the point process of an S4 model object that
defines a finite mixture model. Only available for Poisson mixtures so far.
}
\examples{
\dontrun{
plotPointProc(f_model)
}

}
\seealso{
\code{model}
}
\keyword{internal}
