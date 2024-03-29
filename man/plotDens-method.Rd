% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmcoutputpermhierpost.R
\name{plotDens}
\alias{plotDens}
\title{Plot densities of the parameters and weights}
\arguments{
\item{x}{An \code{mcmcoutput} or \code{mcmcoutputperm} object containing all sampled
values.}

\item{dev}{A logical indicating, if the plots should be shown by a graphical
device. If plots should be stored to a file set \code{dev} to \code{FALSE}.}

\item{...}{Further arguments to be passed to the plotting function.}
}
\value{
Densities of the MCMC samples.
}
\description{
\code{plotDens()} is a class method for \link[=mcmcoutput-class]{mcmcoutput} and
\link[=mcmcoutputperm-class]{mcmcoutputperm} objects. For the former class it
plots densities of MCMC samples and for the latter of the corresponding
permuted samples coming from relabeling.
}
\details{
Calling \code{\link[=plotDens]{plotDens()}} plots densities of the sampled parameters and weights
from MCMC sampling. Note, for relabeled MCMC samples this method is so far
only implemented for mixtures of Poisson and Binomial distributions.
\subsection{Hierarchical priors}{

In case that hierarchical priors had been used in MCMC sampling densities
of the sampled parameters of the hierarchical prior are added to the plot.
}

\subsection{Posterior density parameters}{

In case that posterior density parameters had been stored in MCMC sampling,
densities of these parameters are added to the plot.
}
}
\examples{
# Define a Poisson mixture model with two components.
f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2)
# Simulate data from the mixture model.
f_data <- simulate(f_model)
# Define the hyper-parameters for MCMC sampling.
f_mcmc <- mcmc()
# Define the prior distribution by relying on the data.
f_prior <- priordefine(f_data, f_model)
# Start MCMC sampling.
f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
f_outputperm <- mcmcpermute(f_output)
plotDens(f_outputperm)

}
\seealso{
\itemize{
\item \code{\link[=mixturemcmc]{mixturemcmc()}} for performing MCMC sampling
\item \code{\link[=mcmcpermute]{mcmcpermute()}} for permuting MCMC samples
\item \code{\link[=plotTraces]{plotTraces()}} for plotting the traces of sampled values
\item \code{\link[=plotDens]{plotDens()}} for plotting densities of sampled values
\item \code{\link[=plotSampRep]{plotSampRep()}} for plotting sampling representations of sampled values
\item \code{\link[=plotPointProc]{plotPointProc()}} for plotting point processes for sampled values
\item \code{\link[=plotPostDens]{plotPostDens()}} for plotting the posterior density of component parameters
\item \linkS4class{mcmcoutput} for the class definition of \code{mcmcoutput}
\item \linkS4class{mcmcoutputperm} for the class definition of \code{mcmcoutputperm}
}
}
