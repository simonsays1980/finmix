## Copyright (C) 2013 Lars Simon Zehnder
#
# This file is part of finmix.
#
# finmix is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# finmix is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with finmix. If not, see <http://www.gnu.org/licenses/>.

#' Finmix `mcmcoutputhier` class
#' 
#' @description 
#' This class inherits from the `mcmcoutputbase` class and stores draws from 
#' MCMC sampling with unknown indicators and an hierarchical prior. It adds to 
#' its parent class a slot for storing the parameters of the hierarchical prior.
#' 
#' To use an hierarchical prior in MCMC sampling the `prior` object needs to 
#' have set slot `@hier` to `TRUE`. 
#' 
#' @slot hyper A named list containing the arrays with parameters from the 
#' hierarchical prior.
#' @exportClass mcmcoutputhier
#' @rdname mcmcoutputhier-class
#' 
#' @seealso 
#' * [mcmcoutputbase-class] for the parent class
#' * [prior-class] for the class specifying the prior distribution 
#' * [prior()] for the `prior` class constructor
#' * [priordefine()] for the advanced `prior` class constructor 
.mcmcoutputhier <- setClass("mcmcoutputhier",
  representation(hyper = "list"),
  contains = c("mcmcoutputbase"),
  validity = function(object) {
    ## else: OK
    TRUE
  },
  prototype(hyper = list())
)

#' Shows a summary of an `mcmcoutputhier` object.
#' 
#' @description
#' Calling [show()] on an `mcmcoutputhier` object gives an overview 
#' of the `mcmcoutputhier` object.
#' 
#' @param object An `mcmcoutputhier` object.
#' @returns A console output listing the slots and summary information about
#'   each of them. 
#' @exportMethod show
#' @describeIn mcmcoutput_class
setMethod(
  "show", "mcmcoutputhier",
  function(object) {
    cat("Object 'mcmcoutput'\n")
    cat("     class       :", class(object), "\n")
    cat("     M           :", object@M, "\n")
    cat("     burnin      :", object@burnin, "\n")
    cat("     ranperm     :", object@ranperm, "\n")
    cat(
      "     par         : List of",
      length(object@par), "\n"
    )
    cat(
      "     log         : List of",
      length(object@log), "\n"
    )
    cat(
      "     hyper       : List of",
      length(object@hyper), "\n"
    )
    cat(
      "     ST          :",
      paste(dim(object@ST), collapse = "x"), "\n"
    )
    if (!all(is.na(object@S))) {
      cat(
        "     S           :",
        paste(dim(object@S), collapse = "x"), "\n"
      )
    }
    cat(
      "     NK          :",
      paste(dim(object@NK), collapse = "x"), "\n"
    )
    cat(
      "     clust       :",
      paste(dim(object@clust), collapse = "x"), "\n"
    )
    cat(
      "     model       : Object of class",
      class(object@model), "\n"
    )
    cat(
      "     prior       : Object of class",
      class(object@prior), "\n"
    )
  }
)

#' Plot traces of MCMC sampling
#' 
#' @description 
#' Calling [plotTraces()] plots the MCMC traces of the mixture log-likelihood 
#' , the mixture log-likelihood of the prior distribution, the log-likelihood 
#' of the complete data posterior, or the weights and parameters if `lik` is 
#' set to `1`.s 
#' 
#' If `lik` is set to `0` the parameters of the components and the posterior 
#' parameters are plotted together with `K-1` weights.
#' 
#' @param x An `mcmcoutputhier` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param lik An integer indicating, if the log-likelihood traces should be 
#'   plotted (default). If set to `0` the traces for the parameters 
#'   and weights are plotted instead. 
#' @param col A logical indicating, if the plot should be colored.
#' @param ... Further arguments to be passed to the plotting function.
#' @return A plot of the traces of the MCMC samples.
#' @exportMethod plotTraces
#' 
#' @examples 
#' \dontrun{
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc(storepost = FALSE)
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' plotTraces(f_output, lik = 0)
#' }
#' 
#' @seealso 
#' * [mixturemcmc()] for performing MCMC sampling
#' * [plotHist()] for plotting histograms of sampled values
#' * [plotDens()] for plotting densities of sampled values
#' * [plotSampRep()] for plotting sampling representations of sampled values
#' * [plotPointProc()] for plotting point processes for sampled values
#' * [plotPostDens()] for plotting the posterior density of component parameters
setMethod(
  "plotTraces", signature(
    x = "mcmcoutputhier",
    dev = "ANY",
    lik = "ANY",
    col = "ANY"
  ),
  function(x, dev = TRUE, lik = 1, col = FALSE, ...) {
    dist <- x@model@dist
    if (lik %in% c(0, 1)) {
      if (dist == "poisson" || dist == "cond.poisson") {
        .traces.Poisson.Base.Hier(x, dev)
      } else if (dist == "binomial") {
        .traces.Binomial.Base(x, dev)
      } else if (dist == "exponential") {
        .traces.Exponential.Base(x, dev)
      } else if (dist == "normal") {
        .traces.Normal.Hier(x, dev)
        .traces.Weights.Base(x, dev, col)
      } else if (dist == "student") {
        .traces.Student.Hier(x, dev)
        .traces.Weights.Base(x, dev, col)
      } else if (dist == "normult") {
        .traces.Normult.Hier(x, dev, col)
        .traces.Weights.Base(x, dev, col)
      } else if (dist == "studmult") {
        .traces.Studmult.Hier(x, dev, col)
        .traces.Weights.Base(x, dev, col)
      }
    }
    if (lik %in% c(1, 2)) {
      ## log ##
      .traces.Log.Base(x, dev)
    }
  }
)

#' Plot histograms of the parameters and weights
#' 
#' @description 
#' Calling [plotHist()] plots histograms of the sampled parameters and weights 
#' from MCMC sampling.More specifically, all component parameters, `K-1` of the 
#' weights and the posterior parameters are considered in the histogram plots. 
#' 
#' @param x An `mcmcoutputhier` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Histograms of the MCMC samples.
#' @exportMethod plotHist
#' 
#' @examples
#' \dontrun{ 
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc(storepost = FALSE)
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' plotHist(f_output)
#' }
#' 
#' @seealso 
#' * [mixturemcmc()] for performing MCMC sampling
#' * [plotTraces()] for plotting the traces of sampled values
#' * [plotDens()] for plotting densities of sampled values
#' * [plotSampRep()] for plotting sampling representations of sampled values
#' * [plotPointProc()] for plotting point processes for sampled values
#' * [plotPostDens()] for plotting the posterior density of component parameters
setMethod(
  "plotHist", signature(
    x = "mcmcoutputhier",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    dist <- x@model@dist
    if (dist == "poisson") {
      .hist.Poisson.Base.Hier(x, dev)
    } else if (dist == "binomial") {
      .hist.Binomial.Base(x, dev)
    }
  }
)

#' Plot densities of the parameters and weights
#' 
#' @description 
#' Calling [plotDens()] plots densities of the sampled parameters and weights 
#' from MCMC sampling.More specifically, all component parameters, `K-1` of the 
#' weights and the posterior parameters are considered in the density plots. 
#' 
#' @param x An `mcmcoutputhier` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Densities of the MCMC samples.
#' @exportMethod plotDens
#' 
#' @examples 
#' \dontrun{
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2, 
#'                  indicfix = TRUE)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc(storepost = FALSE)
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' plotDens(f_output)
#' }
#' 
#' @seealso 
#' * [mixturemcmc()] for performing MCMC sampling
#' * [plotTraces()] for plotting the traces of sampled values
#' * [plotHist()] for plotting histograms of sampled values
#' * [plotSampRep()] for plotting sampling representations of sampled values
#' * [plotPointProc()] for plotting point processes for sampled values
#' * [plotPostDens()] for plotting the posterior density of component parameters
setMethod(
  "plotDens", signature(
    x = "mcmcoutputhier",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    dist <- x@model@dist
    if (dist == "poisson") {
      .dens.Poisson.Base.Hier(x, dev)
    } else if (dist == "binomial") {
      .dens.Binomial.Base(x, dev)
    }
  }
)

#' Plot point processes of the component parameters
#' 
#' @description 
#' Calling [plotPointProc()] plots point processes of the sampled component 
#' parameters from MCMC sampling.  
#' 
#' @param x An `mcmcoutputhier` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Point process of the MCMC samples.
#' @exportMethod plotPointProc
#' 
#' @examples 
#' \dontrun{
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc(storepost = FALSE)
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' plotPointProc(f_output)
#' }
#' 
#' @seealso 
#' * [mixturemcmc()] for performing MCMC sampling
#' * [plotTraces()] for plotting the traces of sampled values
#' * [plotHist()] for plotting histograms of sampled values
#' * [plotDens()] for plotting densities of sampled values
#' * [plotSampRep()] for plotting sampling representations of sampled values
#' * [plotPostDens()] for plotting posterior densities for sampled values
setMethod(
  "plotPointProc", signature(
    x = "mcmcoutputhier",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    ## Call 'plotPointProc()' from 'mcmcoutputbase'
    callNextMethod(x, dev, ...)
  }
)

#' Plot sampling representations for the component parameters.
#' 
#' @description 
#' Calling [plotSampRep()] plots sampling representations of the sampled 
#' component parameters from MCMC sampling.  
#' 
#' @param x An `mcmcoutputhier` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Sampling representation of the MCMC samples.
#' @exportMethod plotSampRep
#' 
#' @examples 
#' \dontrun{
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc(storepost = FALSE)
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' plotSampRep(f_output)
#' }
#' 
#' @seealso 
#' * [mixturemcmc()] for performing MCMC sampling
#' * [plotTraces()] for plotting the traces of sampled values
#' * [plotHist()] for plotting histograms of sampled values
#' * [plotDens()] for plotting densities of sampled values
#' * [plotPointProc()] for plotting point processes of sampled values
#' * [plotPostDens()] for plotting posterior densities for sampled values
setMethod(
  "plotSampRep", signature(
    x = "mcmcoutputhier",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    ## Call 'plotSampRep()' from 'mcmcoutputbase'
    callNextMethod(x, dev, ...)
  }
)

#' Plot posterior densities of the component parameters
#' 
#' @description 
#' Calling [plotPostDens()] plots posterior densities of the sampled component 
#' parameters from MCMC sampling, if the number of components is two. 
#' 
#' @param x An `mcmcoutputhier` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Posterior densities of the MCMC samples.
#' @exportMethod plotPostDens
#' @describeIn mcmcoutput_class
#' 
#' @examples 
#' \dontrun{
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc(storepost = FALSE)
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' plotPostDens(f_output)
#' }
#' 
#' @seealso 
#' * [mixturemcmc()] for performing MCMC sampling
#' * [plotTraces()] for plotting the traces of sampled values
#' * [plotHist()] for plotting histograms of sampled values
#' * [plotDens()] for plotting densities of sampled values
#' * [plotSampRep()] for plotting sampling representations of sampled values
#' * [plotPointProc()] for plotting point processes for sampled values
setMethod(
  "plotPostDens", signature(
    x = "mcmcoutputhier",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    ## Call 'plotPostDens()' from 'mcmcoutputbase'
    callNextMethod(x, dev, ...)
  }
)

#' Constructs a sub-chain of MCMC samples 
#' 
#' @description 
#' Calling [subseq()] constructs an MCMC sub-chain from the samples in the 
#' passed-in `mcmcoutput` object specfied by the index `array` in `index`. This 
#' can be advantageous, if chains are non-stationary. For successful MCMC 
#' sampling the chain must be converged to the target distribution, the true 
#' distribution of parameters, weights and indicators.
#' 
#' Note, this method calls the equivalent method of the parent class and then 
#' adds to it the sub-chains for the parameters of the hierarchical prior.
#' 
#' @param object An `mcmcoutput` object containing all sampled values.
#' @param index An array specifying the extraction of the sub-chain.
#' @return An `mcmcoutput` object containing the values from the sub-chain.
#' @noRd
setMethod(
  "subseq", signature(
    object = "mcmcoutputhier",
    index = "array"
  ),
  function(object, index) {
    ## Call 'subseq()' method from 'mcmcoutputfixhier'
    as(object, "mcmcoutputbase") <- callNextMethod(object, index)
    dist <- object@model@dist
    if (dist == "poisson") {
      .subseq.Poisson.Hier(object, index)
    } else if (dist %in% c("normal", "student")) {
      .subseq.Norstud.Hier(object, index)
    } else if (dist %in% c("normult", "studmult")) {
      .subseq.Normultstud.Hier(object, index)
    }
  }
)

#' Swaps elements between components
#' 
#' @description 
#' Not yet implemented.
#' 
#' @param object An `mcmcoutput` object containing the sampled values.
#' @param index An array specifying the extraction of the values.
#' @return An `mcmcoutput` object with swapped elements.
#' @noRd 
setMethod(
  "swapElements", signature(
    object = "mcmcoutputhier",
    index = "array"
  ),
  function(object, index) {
    ## Check arguments, TODO: .validObject ##
    ## Call method 'swapElements()' from 'mcmcoutputbase'
    callNextMethod(object, index)
  }
)

#' Getter method of `mcmcoutputhier` class.
#' 
#' Returns the `hyper` slot.
#' 
#' @param object An `mcmcoutputhier` object.
#' @returns The `hyper` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc(storepost = FALSE)
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' # Get the slot.
#' getHyper(f_output)
#' 
#' @seealso 
#' * [mcmcoutputhier-class] for the class definition
#' * [mixturemcmc()] for performing MCMC sampling
setMethod(
  "getHyper", "mcmcoutputhier",
  function(object) {
    return(object@hyper)
  }
)


## No setters as users are not intended to manipulate ##
## this object ##

### Private functions.
### These functions are not exported.

### Plot
### Plot traces

#' Plots traces of Poisson mixture samples
#' 
#' @description 
#' For internal usage only. This function plots the traces for sampled values 
#' from a Poisson mixture model.
#' 
#' @param x An `mcmcoutputhier` object containing all samples.
#' @param dev A logical indicating if the plot should be shown by a grapical 
#'   device.
#' @return A plot of the traces of sampled values.
#' @noRd
#' 
#' @seealso 
#' * [plotTraces()] for the calling function
".traces.Poisson.Base.Hier" <- function(x, dev) {
  K <- x@model@K
  trace.n <- K * 2
  if (.check.grDevice() && dev) {
    dev.new(title = "Traceplots")
  }
  par(
    mfrow = c(trace.n, 1), mar = c(1, 0, 0, 0),
    oma = c(4, 5, 4, 4)
  )
  lambda <- x@par$lambda
  for (k in 1:K) {
    plot(lambda[, k],
      type = "l", axes = F,
      col = "gray20", xlab = "", ylab = ""
    )
    axis(2, las = 2, cex.axis = 0.7)
    mtext(
      side = 2, las = 2, bquote(lambda[k = .(k)]),
      cex = 0.6, line = 3
    )
  }
  weight <- x@weight
  for (k in 1:(K - 1)) {
    plot(weight[, k],
      type = "l", axes = F,
      col = "gray47", xlab = "", ylab = ""
    )
    axis(2, las = 2, cex.axis = 0.7)
    mtext(
      side = 2, las = 2, bquote(eta[k = .(k)]),
      cex = 0.6, line = 3
    )
  }
  b <- x@hyper$b
  plot(b,
    type = "l", axes = F,
    col = "gray68", xlab = "", ylab = ""
  )
  axis(2, las = 2, cex.axis = 0.7)
  mtext(side = 2, las = 2, "b", cex = 0.6, line = 3)
  axis(1)
  mtext(side = 1, "Iterations", cex = 0.7, line = 3)
}

### Plot Histograms

#' Plot histograms of Poisson samples
#' 
#' @description 
#' For internal usage only. This function plots histograms of sampled Poisson 
#' parameters and weights. In addition it plots the histogram of the 
#' parameter `b` of the hierarchical prior.
#' 
#' @param x An `mcmcoutputhier` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with histograms for the smapled parameters and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotHist()] for the calling function
".hist.Poisson.Base.Hier" <- function(x, dev) {
  K <- x@model@K
  if (.check.grDevice() && dev) {
    dev.new(title = "Histograms")
  }
  lambda <- x@par$lambda
  weight <- x@weight
  b <- x@hyper$b
  vars <- cbind(lambda, weight[, seq(1, K - 1)], b)
  lab.names <- vector("list", 2 * K)
  for (k in seq(1, K)) {
    lab.names[[k]] <- bquote(lambda[.(k)])
  }
  for (k in seq(K + 1, 2 * K - 1)) {
    lab.names[[k]] <- bquote(eta[.(k - K)])
  }
  lab.names[[2 * K]] <- "b"
  .symmetric.Hist(vars, lab.names)
}

### Plot Densities

#' Plot densities of Poisson samples
#' 
#' @description 
#' For internal usage only. This function plots densities of sampled Poisson 
#' parameters and weights. In addition it plots the Kernel densities of the 
#' parameter `b` of the hierarchical prior.
#' 
#' @param x An `mcmcoutputhier` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with densities for the sampled parameters and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotDens()] for the calling function
".dens.Poisson.Base.Hier" <- function(x, dev) {
  K <- x@model@K
  if (.check.grDevice() && dev) {
    dev.new(title = "Densities")
  }
  lambda <- x@par$lambda
  weight <- x@weight
  b <- x@hyper$b
  vars <- cbind(lambda, weight[, seq(1, K - 1)], b)
  lab.names <- vector("list", 2 * K)
  for (k in seq(1, K)) {
    lab.names[[k]] <- bquote(lambda[.(k)])
  }
  for (k in seq(K + 1, 2 * K - 1)) {
    lab.names[[k]] <- bquote(eta[.(k - K)])
  }
  lab.names[[2 * K]] <- "b"
  .symmetric.Dens(vars, lab.names)
}
