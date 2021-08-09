\name{simulate}
\alias{simulate}
\title{
Simulate data from a finite mixture model
}
\description{
Simulate data from a finite mixture model defined by a finmix 'model' object. Simulated data includes values and
indicators (allocations). The return value is an object of class \code{fdata} defined in the \code{finmix} package.
}
\usage{
simulate(model, N = 100, varargin, seed = 0)
}

\arguments{
  \item{model}{
    an S4 object of class \code{model} defined in the \code{finmix} package with fully specified parameters in 
    \code{@par}, weights in \code{@weight}, number of components in \code{@K} and distribution in
    \code{@dist}.
}
  \item{N}{
    number of values to be simulated. If \code{typeof(n) == "numeric"} \code{n} is rounded to the next 
    integer value.
}
  \item{varargin}{
    an S4 object of class \code{fdata} defined in the \code{finmix} package; optional.
  }
  \item{seed}{
    The seed for the random generator.
}
}
\details{
    This function is an \code{S4} method of the \code{model} class in \code{finmix}. All slots of \code{model} 
    must be specified. 
    The function returns an \code{fdata} object with indicators in slot \code{@S} and values in slot 
    \code{@y}. Slot \code{@sim} is set to \code{TRUE} to indicate that the data was simulated. If an 
    \code{fdata} object is already provided in the argument \code{varargin}, this object gets the simulated 
    data assigned and is returned. 
    For a Binomial mixture model, i.e. \code{@dist == "binomial"} in \code{model}, repetitions are needed. 
    If not provided by \code{varargin} in slot \code{@T} repetitions are set all to 1. 
}
\value{
    An S4 object of class \code{fdata} defined by the \code{finmix} package.
}

\note{ 
    This function is a class method of the S4 \code{model} class defined in the \code{finmix} package.
}


\seealso{
\code{\link{model}},  \code{\link{fdata}},  \code{\link{finmix}}
}
\examples{
\dontrun{

## Simulate a Poisson model with 2 components
fm.model <- model(dist = "poisson", K = 2)
weight <- c(.3, .7)
setWeight(fm.model) <- weight
pars <- list(lambda = c(312, 80))
setPar(fm.model) <- pars
fm.fdata <- simulate(fm.model, N = 1000, seed = 123456)

## Simulate a Binomial model with 3 components
fm.model <- model(dist = "binomial", K = 3,
                  weight = c(.4, .6),
                  par = list(p = c(.2, .9)))
fm.fdata <- fdata(T = matrix(seq(1, 1000)))
RNGkind("Marsaglia-Multicarry")
sim.fdata <- simulate(fm.model, N = 1000, fm.fdata)
   
}
}
