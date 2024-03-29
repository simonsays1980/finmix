% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmcoutputpermhierpost.R
\name{swapElements}
\alias{swapElements}
\title{Swap elements of MCMC samples}
\arguments{
\item{object}{An \code{mcmcoutput} object containing the
sampled values.}

\item{index}{An array specifying the extraction of the values.}
}
\value{
An \code{mcmcoutput} object with swapped elements.
}
\description{
Calling \code{swapElements()} on an \code{mcmcoutput} object
swaps all labels by the schema given in the \code{index} argument.
}
\details{
This function is merely a utility function that simplifies relabeling for
users and developers. For relabeling the labels have to be permuted and
depending on the MCMC sampling chosen there could be a lot of different
slots that need to be permuted. \code{swapElements()} swaps the elements in any
slot that needs to be relabeled.
}
\examples{
\dontrun{
# Generate a model of Poisson distributions.
f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2)
# Simulate data from the model.
f_data <- simulate(f_model)
# Define the hyper-parameters for MCMC sampling.
f_mcmc <- mcmc(storepost = FALSE)
# Define the prior distribution by relying on the data.
f_prior <- priordefine(f_data, f_model)
# Start MCMC sampling.
f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
index <- matrix(c(1, 2), nrow = getM(f_output) + 1, 
                ncol = 2)[1:getM(f_output),]
swapElements(f_output, index)
}

}
\seealso{
\itemize{
\item \linkS4class{mcmcoutput} for the class definition
\item \code{\link[=subseq]{subseq()}} for generating sub-chains from MCMC samples
\item \code{\link[=mcmcpermute]{mcmcpermute()}} for a calling function
}
}
