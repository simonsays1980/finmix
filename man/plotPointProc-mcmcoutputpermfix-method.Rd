% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmcoutputpermfix.R
\name{plotPointProc,mcmcoutputpermfix-method}
\alias{plotPointProc,mcmcoutputpermfix-method}
\title{Plot point processes of the component parameters}
\usage{
\S4method{plotPointProc}{mcmcoutputpermfix}(x, dev = TRUE, ...)
}
\arguments{
\item{x}{An \code{mcmcoutputpermfix} object containing all sampled values.}

\item{dev}{A logical indicating, if the plots should be shown by a graphical
device. If plots should be stored to a file set \code{dev} to \code{FALSE}.}

\item{...}{Further arguments to be passed to the plotting function.}
}
\value{
Densities of the MCMC samples.
}
\description{
Calling \code{\link[=plotPointProc]{plotPointProc()}} plots point processes of the sampled component
parameters from MCMC sampling.
}
\examples{
\dontrun{
# Define a Poisson mixture model with two components.
f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2, 
                 indicfix = TRUE)
# Simulate data from the mixture model.
f_data <- simulate(f_model)
# Define the hyper-parameters for MCMC sampling.
f_mcmc <- mcmc(storepost = FALSE)
# Define the prior distribution by relying on the data.
f_prior <- priordefine(f_data, f_model)
# Do not use a hierarchical prior.
setHier(f_prior) <- FALSE
# Start MCMC sampling.
f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
f_outputperm <- mcmcpermute(f_output)
plotPointProc(f_outputperm)
}

}
\seealso{
\itemize{
\item \code{\link[=mixturemcmc]{mixturemcmc()}} for performing MCMC sampling
\item \code{\link[=mcmcpermute]{mcmcpermute()}} for permuting MCMC samples
\item \code{\link[=plotTraces]{plotTraces()}} for plotting the traces of sampled values
\item \code{\link[=plotHist]{plotHist()}} for plotting histograms of sampled values
\item \code{\link[=plotDens]{plotDens()}} for plotting densities of sampled values
\item \code{\link[=plotSampRep]{plotSampRep()}} for plotting sampling representations of sampled values
\item \code{\link[=plotPostDens]{plotPostDens()}} for plotting posterior densities for sampled values
}
}
\keyword{internal}
