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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with finmix. If not, see <http://www.gnu.org/licenses/>.

# TODO: CHange examples to storepost = FALSE

#' Finmix `mcmcoutput` base class for unknown indicators
#' 
#' @description
#' This class defines the basic slots for the MCMC sampling output when 
#' indicators are not known. It inherits from the 
#' [mcmcoutputfix-class]. 
#' 
#' @slot weight An `array` of dimension `M x K` containing the sampled 
#'   weights.
#' @slot entropy An `array` of dimension `M x 1` containing the entropy 
#'   for each MCMC draw.
#' @slot ST An `array` of dimension `M x 1` containing all MCMC states, 
#'   for the last observation in slot `y` of the `fdata` object passed in to 
#'   [mixturemcmc()] where a state is defined for non-Markov models as the 
#'   last indicator of this observation. 
#' @slot S An `array` of dimension `N x storeS` containing the last 
#'   `storeS` indicators sampled. `storeS` is defined in the slot `storeS` of 
#'   the `mcmc` object passed into [mixturemcmc()].
#' @slot NK An `array` of dimension `M x K` containing the number of 
#'   observations assigned to each component for each MCMC draw.
#' @slot clust An `array` of dimension `N x 1` containing the recent 
#'   indicators defining the last "clustering" of observations into the 
#'   mixture components.
#' @exportClass mcmcoutputbase
#' @rdname mcmcoutputbase-class
.mcmcoutputbase <- setClass("mcmcoutputbase",
  representation(
    weight = "array",
    entropy = "array",
    ST = "array",
    S = "array",
    NK = "array",
    clust = "array"
  ),
  contains = c("mcmcoutputfix"),
  validity = function(object) {
    ## else: OK
    TRUE
  },
  prototype(
    weight = array(),
    entropy = array(),
    ST = array(),
    S = array(),
    NK = array(),
    clust = array()
  )
)

#' Shows a summary of an `mcmcoutputbase` object.
#' 
#' Calling [show()] on an `mcmcoutputbase` object gives an overview 
#' of the `mcmcoutputbase` object.
#' 
#' @param object An `mcmcoutputbase` object.
#' @returns A console output listing the slots and summary information about
#'   each of them. 
#' @exportMethod show
#' @noRd 
setMethod(
  "show", "mcmcoutputbase",
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
      "     entropy     : ",
        paste(dim(object@entropy), collapse = "x"), "\n"
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
#' set to `0`. 
#' 
#' If `lik` is set to `0` the parameters of the components and the posterior 
#' parameters are plotted together with `K-1` weights.
#' 
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param lik An integer indicating, if the log-likelihood traces should be 
#'   plotted (default). If set to `0` the traces for the parameters 
#'   and weights are plotted instead. 
#' @param col A logical indicating, if the plot should be colored.
#' @param ... Further arguments to be passed to the plotting function.
#' @return A plot of the traces of the MCMC samples.
#' @exportMethod plotTraces
#' @noRd
#' 
#' @examples 
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc()
#' # Complete object slots for consistency. 
#' (f_data ~ f_model ~ f_mcmc) %=% mcmcstart(f_data, f_model, f_mcmc)
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' plotTraces(f_output, lik = 0)
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
    x = "mcmcoutputbase",
    dev = "ANY",
    lik = "ANY",
    col = "ANY"
  ),
  function(x, dev = TRUE, lik = 1, col = FALSE, ...) {
    dist <- x@model@dist
    if (lik %in% c(0, 1)) {
      if (dist == "poisson" || dist == "cond.poisson") {
        .traces.Poisson.Base(x, dev)
      } else if (dist == "binomial") {
        .traces.Binomial.Base(x, dev)
      } else if (dist == "exponential") {
        .traces.Exponential.Base(x, dev)
      } else if (dist == "normal") {
        .traces.Normal(x, dev)
        .traces.Weights.Base(x, dev, col)
      } else if (dist == "student") {
        .traces.Student(x, dev)
        .traces.Weights.Base(x, dev, col)
      } else if (dist == "normult") {
        .traces.Normult(x, dev, col)
        .traces.Weights.Base(x, dev, col)
      } else if (dist == "studmult") {
        .traces.Studmult(x, dev, col)
        .traces.Weights.Base(x, dev, col)
      }
    }
    if (lik %in% c(1, 2)) {
      ## log ##
      .traces.Log.Base(x, dev, col)
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
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Histograms of the MCMC samples.
#' @exportMethod plotHist
#' @noRd
#' 
#' @examples 
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc()
#' # Complete object slots for consistency. 
#' (f_data ~ f_model ~ f_mcmc) %=% mcmcstart(f_data, f_model, f_mcmc)
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' plotHist(f_output)
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
    x = "mcmcoutputbase",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    dist <- x@model@dist
    if (dist == "poisson") {
      .hist.Poisson.Base(x, dev)
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
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Densities of the MCMC samples.
#' @exportMethod plotDens
#' @noRd
#' 
#' @examples 
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc()
#' # Complete object slots for consistency. 
#' (f_data ~ f_model ~ f_mcmc) %=% mcmcstart(f_data, f_model, f_mcmc)
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' plotDens(f_output)
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
    x = "mcmcoutputbase",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    dist <- x@model@dist
    if (dist == "poisson") {
      .dens.Poisson.Base(x, dev)
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
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Point processes of the MCMC samples.
#' @exportMethod plotPointProc
#' @noRd
#' 
#' @examples 
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc()
#' # Complete object slots for consistency. 
#' (f_data ~ f_model ~ f_mcmc) %=% mcmcstart(f_data, f_model, f_mcmc)
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' plotPointProc(f_output)
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
    x = "mcmcoutputbase",
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
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Sampling representations of the MCMC samples.
#' @exportMethod plotSampRep
#' @noRd
#' 
#' @examples 
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc()
#' # Complete object slots for consistency. 
#' (f_data ~ f_model ~ f_mcmc) %=% mcmcstart(f_data, f_model, f_mcmc)
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' plotSampRep(f_output)
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
    x = "mcmcoutputbase",
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
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Posterior densities of the MCMC samples.
#' @exportMethod plotPostDens
#' @noRd
#' 
#' @examples 
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc()
#' # Complete object slots for consistency. 
#' (f_data ~ f_model ~ f_mcmc) %=% mcmcstart(f_data, f_model, f_mcmc)
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' plotPostDens(f_output)
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
    x = "mcmcoutputbase",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    ## Call 'plotPostDens' from 'mcmcoutputfixhier'
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
#' @param object An `mcmcoutput` object containing all sampled values.
#' @param index An array specifying the extraction of the sub-chain.
#' @return An `mcmcoutput` object containing the values from the sub-chain.
#' @noRd
setMethod(
  "subseq", signature(
    object = "mcmcoutputbase",
    index = "array"
  ),
  function(object, index) {
    ## Call 'subseq()' method from 'mcmcoutputfix'
    as(object, "mcmcoutputfix") <- callNextMethod(object, index)
    ## Change owned slots ##
    .subseq.Base(object, index)
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
#' @keywords internal 
setMethod(
  "swapElements", signature(
    object = "mcmcoutputbase",
    index = "array"
  ),
  function(object, index) {
    if (object@model@K == 1) {
      return(object)
    } else {
      ## Call method 'swapElements()' from 'mcmcoutputfix'
      as(object, "mcmcoutputfix") <- callNextMethod(object, index)
      .swapElements.Base(object, index)
    }
  }
)

#' Getter method of `mcmcoutput` class.
#' 
#' Returns the `weight` slot.
#' 
#' @param object An `mcmcoutput` object.
#' @returns The `weight` slot of the `object`.
#' @exportMethod getWeight
#' @noRd
#' 
#' @examples 
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc()
#' # Complete object slots for consistency. 
#' (f_data ~ f_model ~ f_mcmc) %=% mcmcstart(f_data, f_model, f_mcmc)
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' # Get the slot.
#' getWeight(f_output)
#' 
#' @seealso 
#' * [mcmcoutput][mcmcoutputbase-class] for the class definition
#' * [mixturemcmc()] for performing MCMC sampling
setMethod(
  "getWeight", "mcmcoutputbase",
  function(object) {
    return(object@weight)
  }
)

#' Getter method of `mcmcoutput` class.
#' 
#' Returns the `entropy` slot.
#' 
#' @param object An `mcmcoutput` object.
#' @returns The `entropy` slot of the `object`.
#' @exportMethod getEntropy
#' @noRd
#' 
#' @examples 
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc()
#' # Complete object slots for consistency. 
#' (f_data ~ f_model ~ f_mcmc) %=% mcmcstart(f_data, f_model, f_mcmc)
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' # Get the slot.
#' getEntropy(f_output)
#' 
#' @seealso 
#' * [mcmcoutput][mcmcoutputbase-class] for the class definition
#' * [mixturemcmc()] for performing MCMC sampling
setMethod(
  "getEntropy", "mcmcoutputbase",
  function(object) {
    return(object@entropy)
  }
)

#' Getter method of `mcmcoutput` class.
#' 
#' Returns the `ST` slot.
#' 
#' @param object An `mcmcoutput` object.
#' @returns The `ST` slot of the `object`.
#' @exportMethod getST
#' @noRd
#' 
#' @examples 
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc()
#' # Complete object slots for consistency. 
#' (f_data ~ f_model ~ f_mcmc) %=% mcmcstart(f_data, f_model, f_mcmc)
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' # Get the slot.
#' getST(f_output)
#' 
#' @seealso 
#' * [mcmcoutput][mcmcoutputbase-class] for the class definition
#' * [mixturemcmc()] for performing MCMC sampling
setMethod(
  "getST", "mcmcoutputbase",
  function(object) {
    return(object@ST)
  }
)

#' Getter method of `mcmcoutput` class.
#' 
#' Returns the `S` slot.
#' 
#' @param object An `mcmcoutput` object.
#' @returns The `S` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc()
#' # Complete object slots for consistency. 
#' (f_data ~ f_model ~ f_mcmc) %=% mcmcstart(f_data, f_model, f_mcmc)
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' # Get the slot.
#' getS(f_output)
#' 
#' @seealso 
#' * [mcmcoutput][mcmcoutputbase-class] for the class definition
#' * [mixturemcmc()] for performing MCMC sampling
setMethod(
  "getS", "mcmcoutputbase",
  function(object) {
    return(object@S)
  }
)

#' Getter method of `mcmcoutput` class.
#' 
#' Returns the `NK` slot.
#' 
#' @param object An `mcmcoutput` object.
#' @returns The `NK` slot of the `object`.
#' @exportMethod getNK
#' @noRd
#' 
#' @examples 
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc()
#' # Complete object slots for consistency. 
#' (f_data ~ f_model ~ f_mcmc) %=% mcmcstart(f_data, f_model, f_mcmc)
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' # Get the slot.
#' getNK(f_output)
#' 
#' @seealso 
#' * [mcmcoutput][mcmcoutputbase-class] for the class definition
#' * [mixturemcmc()] for performing MCMC sampling
setMethod(
  "getNK", "mcmcoutputbase",
  function(object) {
    return(object@NK)
  }
)

#' Getter method of `mcmcoutput` class.
#' 
#' Returns the `clust` slot.
#' 
#' @param object An `mcmcoutput` object.
#' @returns The `clust` slot of the `object`.
#' @exportMethod getClust
#' @noRd
#' 
#' @examples 
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc()
#' # Complete object slots for consistency. 
#' (f_data ~ f_model ~ f_mcmc) %=% mcmcstart(f_data, f_model, f_mcmc)
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' # Get the slot.
#' getClust(f_output)
#' 
#' @seealso 
#' * [mcmcoutput][mcmcoutputbase-class] for the class definition
#' * [mixturemcmc()] for performing MCMC sampling
setMethod(
  "getClust", "mcmcoutputbase",
  function(object) {
    return(object@clust)
  }
)

## No setters as users are not intended to manipulate ##
## this object. ##

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
#' @param x An `mcmcoutput` object containing all samples.
#' @param dev A logical indicating if the plot should be shown by a grapical 
#'   device.
#' @return A plot of the traces of sampled values.
#' @noRd
#' 
#' @seealso 
#' * [plotTraces()] for the calling function
".traces.Poisson.Base" <- function(x, dev) {
  K <- x@model@K
  trace.n <- K * 2 - 1
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
    axis(2, las = 2, cex.axis = .7)
    mtext(
      side = 2, las = 2, bquote(lambda[k = .(k)]),
      cex = .6, line = 3
    )
  }
  weight <- x@weight
  for (k in 1:(K - 1)) {
    plot(weight[, k],
      type = "l", axes = F,
      col = "gray47", xlab = "", ylab = ""
    )
    axis(2, las = 2, cex.axis = .7)
    mtext(
      side = 2, las = 2, bquote(eta[k = .(k)]),
      cex = .6, line = 3
    )
  }
  axis(1)
  mtext(side = 1, "Iterations", cex = .7, line = 3)
}

#' Plots traces of Binomial mixture samples
#' 
#' @description 
#' For internal usage only. This function plots the traces for sampled values 
#' from a Binomial mixture model.
#' 
#' @param x An `mcmcoutput` object containing all samples.
#' @param dev A logical indicating if the plot should be shown by a grapical 
#'   device.
#' @return A plot of the traces of sampled values.
#' @noRd
#' 
#' @seealso 
#' * [plotTraces()] for the calling function
".traces.Binomial.Base" <- function(x, dev) {
  K <- x@model@K
  trace.n <- K * 2 - 1
  if (.check.grDevice() && dev) {
    dev.new(title = "Traceplots")
  }
  par(
    mfrow = c(trace.n, 1), mar = c(1, 0, 0, 0),
    oma = c(4, 5, 4, 4)
  )
  p <- x@par$p
  for (k in 1:K) {
    plot(p[, k],
      type = "l", axes = F, col = "gray20",
      xlab = "", ylab = ""
    )
    axis(2, las = 2, cex.axis = .7)
    mtext(
      side = 2, las = 2, bquote(p[k = .(k)]),
      cex = .6, line = 3
    )
  }
  weight <- x@weight
  for (k in 1:(K - 1)) {
    plot(weight[, k],
      type = "l", axes = F,
      col = "gray47", xlab = "", ylab = ""
    )
    axis(2, las = 2, cex.axis = .7)
    mtext(
      side = 2, las = 2, bquote(eta[k = .(k)]),
      cex = .6, line = 3
    )
  }
  axis(1)
  mtext(side = 1, "Iterations", cex = .7, line = 3)
}

#' Plots traces of exponential mixture samples
#' 
#' @description 
#' For internal usage only. This function plots the traces for sampled values 
#' from a exponential mixture model.
#' 
#' @param x An `mcmcoutput` object containing all samples.
#' @param dev A logical indicating if the plot should be shown by a grapical 
#'   device.
#' @return A plot of the traces of sampled values.
#' @noRd
#' 
#' @seealso 
#' * [plotTraces()] for the calling function
".traces.Exponential.Base" <- function(x, dev) {
  K <- x@model@K
  trace.n <- K * 2 - 1
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
    axis(2, las = 2, cex.axis = .7)
    mtext(
      side = 2, las = 2, bquote(lambda[k = .(k)]),
      cex = .6, line = 3
    )
  }
  weight <- x@weight
  for (k in 1:(K - 1)) {
    plot(weight[, k],
      type = "l", axes = F,
      col = "gray47", xlab = "", ylab = ""
    )
    axis(2, las = 2, cex.axis = .7)
    mtext(
      side = 2, las = 2, bquote(eta[k = .(k)]),
      cex = .6, line = 3
    )
  }
  axis(1)
  mtext(side = 1, "Iterations", cex = .7, line = 3)
}

#' Plots traces of weights from any mixture model
#' 
#' @description 
#' For internal usage only. This function plots the traces for sampled values 
#' from any mixture model.
#' 
#' @param x An `mcmcoutput` object containing all samples.
#' @param dev A logical indicating if the plot should be shown by a grapical 
#'   device.
#' @return A plot of the traces of sampled values.
#' @noRd
#' 
#' @seealso 
#' * [plotTraces()] for the calling function
".traces.Weights.Base" <- function(x, dev, col) {
  weight <- x@weight
  K <- x@model@K
  if (.check.grDevice() && dev) {
    dev.new(title = "Traceplots Weights")
  }
  if (col) {
    cscale <- rainbow(K, start = 0.5, end = 0)
  } else {
    cscale <- gray.colors(K, start = 0.5, end = 0.15)
  }

  plot(weight[, 1],
    type = "l", axes = F,
    col = cscale[1], xlab = "", ylab = "", ylim = c(0, 1.2)
  )
  for (k in 2:K) {
    lines(weight[, k], col = cscale[k])
  }
  axis(2, las = 2, cex.axis = .7)
  mtext(
    side = 2, las = 2, bquote(eta),
    cex = .6, line = 3
  )
  name <- vector("character", K)
  for (k in 1:K) {
    name[k] <- paste("k = ", k, sep = "")
  }
  legend("top",
    legend = name, col = cscale, lty = 1,
    horiz = TRUE, cex = .7
  )
  axis(1)
  mtext(side = 1, "Iterations", cex = .7, line = 3)
}

#' Plots traces of log-likelihood samples
#' 
#' @description 
#' For internal usage only. This function plots the traces for sampled values 
#' from any mixture model.
#' 
#' @param x An `mcmcoutput` object containing all samples.
#' @param dev A logical indicating if the plot should be shown by a grapical 
#'   device.
#' @return A plot of the traces of sampled values.
#' @noRd
#' 
#' @seealso 
#' * [plotTraces()] for the calling function
".traces.Log.Base" <- function(x, dev, col = FALSE) {
  if (.check.grDevice() && dev) {
    dev.new(title = "Log Likelihood Traceplots")
  }
  if (col) {
    cscale <- rainbow(3, start = 0.5, end = 0)
  } else {
    cscale <- gray.colors(3, start = 0.5, end = 0.15)
  }
  par(
    mfrow = c(3, 1), mar = c(1, 0, 0, 0),
    oma = c(4, 5, 4, 4)
  )
  mixlik <- x@log$mixlik
  plot(mixlik,
    type = "l", axes = F,
    col = cscale[3], xlab = "", ylab = ""
  )
  axis(2, las = 2, cex.axis = 0.7)
  mtext(
    side = 2, las = 3, "mixlik", cex = 0.6,
    line = 3
  )
  mixprior <- x@log$mixprior
  plot(mixprior,
    type = "l", axes = F,
    col = cscale[2], xlab = "", ylab = ""
  )
  axis(2, las = 2, cex.axis = 0.7)
  mtext(
    side = 2, las = 3, "mixprior", cex = 0.6,
    line = 3
  )
  cdpost <- x@log$cdpost
  plot(cdpost,
    type = "l", axes = F,
    col = cscale[1], xlab = "", ylab = ""
  )
  axis(2, las = 2, cex.axis = 0.7)
  mtext(
    side = 2, las = 3, "cdpost", cex = 0.6,
    line = 3
  )
  axis(1)
  mtext(side = 1, "Iterations", cex = 0.7, line = 3)
}

### Histograms

#' Plot histograms of Poisson samples
#' 
#' @description 
#' For internal usage only. This function plots histograms of sampled Poisson 
#' parameters and weights.
#' 
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with histograms for the smapled parameters and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotHist()] for the calling function
".hist.Poisson.Base" <- function(x, dev) {
  K <- x@model@K
  if (.check.grDevice() && dev) {
    dev.new(title = "Histograms")
  }
  lambda <- x@par$lambda
  weight <- x@weight
  vars <- cbind(lambda, weight[, seq(1, K - 1)])
  lab.names <- vector("list", 2 * K - 1)
  for (k in seq(1, K)) {
    lab.names[[k]] <- bquote(lambda[.(k)])
  }
  for (k in seq(K + 1, 2 * K - 1)) {
    lab.names[[k]] <- bquote(eta[.(k - K)])
  }
  .symmetric.Hist(vars, lab.names)
}

#' Plot histograms of Binomial samples
#' 
#' @description 
#' For internal usage only. This function plots histograms of sampled Binomial 
#' parameters and weights.
#' 
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with histograms for the sampled parameters and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotHist()] for the calling function
".hist.Binomial.Base" <- function(x, dev) {
  K <- x@model@K
  if (.check.grDevice() && dev) {
    dev.new(title = "Histograms")
  }
  p <- x@par$p
  weight <- x@weight
  vars <- cbind(p, weight[, seq(1, K - 1)])
  lab.names <- vector("list", 2 * K - 1)
  for (k in seq(1, K)) {
    lab.names[[k]] <- bquote(p[.(k)])
  }
  for (k in seq(K + 1, 2 * K - 1)) {
    lab.names[[k]] <- bquote(eta[.(k - K)])
  }
  .symmetric.Hist(vars, lab.names)
}

#' Plot histograms of exponential samples
#' 
#' @description 
#' For internal usage only. This function plots histograms of sampled 
#' exponential parameters and weights.
#' 
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with histograms for the smapled parameters and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotHist()] for the calling function
".hist.Exponential.Base" <- function(x, dev) {
  K <- x@model@K
  if (.check.grDevice() && dev) {
    dev.new(title = "Histograms")
  }
  lambda <- x@par$lambda
  weight <- x@weight
  vars <- cbind(lambda, weight[, seq(1, K - 1)])
  lab.names <- vector("list", 2 * K - 1)
  for (k in seq(1, K)) {
    lab.names[[k]] <- bquote(lambda[.(k)])
  }
  for (k in seq(K + 1, 2 * K - 1)) {
    lab.names[[k]] <- bquote(eta[.(k - K)])
  }
  .symmetric.Hist(vars, lab.names)
}

#' Plot histograms of normal samples
#' 
#' @description 
#' For internal usage only. This function plots histograms of sampled normal
#' parameters and weights.
#' 
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with histograms for the smapled parameters and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotHist()] for the calling function
".hist.Normal.Base" <- function(x, dev) {
  K <- x@model@K
  .hist.Normal(x, dev)
  if (K > 1) {
    weight <- x@weight
    weights.lab.names <- vector("list", K)
    for (k in 1:K) {
      weights.lab.names[[k]] <- bquote(eta[.(k)])
    }
    if (K > 1) {
      if (.check.grDevice() && dev) {
        dev.new(title = "Histograms Weights")
      }
      .symmetric.Hist(weight, weights.lab.names)
    }
  }
}

#' Plot histograms of Student-t samples
#' 
#' @description 
#' For internal usage only. This function plots histograms of sampled Student-t 
#' parameters and weights.
#' 
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with histograms for the smapled parameters and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotHist()] for the calling function
".hist.Student.Base" <- function(x, dev) {
  K <- x@model@K
  .hist.Student(x, dev)
  if (K > 1) {
    weight <- x@weight
    weight.lab.names <- vector("list", K)
    for (k in 1:K) {
      weight.lab.names[[k]] <- bquote(eta[.(k)])
    }
    if (.check.grDevice() && dev) {
      dev.new(title = "Histograms Weights")
    }
    .symmetric.Hist(weight, weight.lab.names)
  }
}

#' Plot histograms of multivariate normal samples
#' 
#' @description 
#' For internal usage only. This function plots histograms of sampled 
#' multivariate normal parameters and weights.
#' 
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with histograms for the smapled parameters and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotHist()] for the calling function
".hist.Normult.Base" <- function(x, dev) {
  K <- x@model@K
  .hist.Normult(x, dev)
  if (K > 1) {
    weight <- x@weight
    weight.lab.names <- vector("list", K)
    for (k in 1:K) {
      weight.lab.names[[k]] <- bquote(eta[.(k)])
    }
    if (.check.grDevice() && dev) {
      dev.new(title = "Histograms Weights")
    }
    .symmetric.Hist(weight, weight.lab.names)
  }
}

#' Plot histograms of multivariate Student-t samples
#' 
#' @description 
#' For internal usage only. This function plots histograms of sampled 
#' multivariate Student-t parameters and weights.
#' 
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with histograms for the smapled parameters and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotHist()] for the calling function
".hist.Studmult.Base" <- function(x, dev) {
  K <- x@model@K
  .hist.Studmult(x, dev)
  if (K > 1) {
    weight <- x@weight
    weight.lab.names <- vector("list", K)
    for (k in 1:K) {
      weight.lab.names[[k]] <- bquote(eta[.(k)])
    }
    .symmetric.Hist(weight, weight.lab.names)
  }
}

### Densities
### Densities Poisson: Plots Kernel densities for the Poisson
### parameters and the weights.
#' Plot densities of Poisson samples
#' 
#' @description 
#' For internal usage only. This function plots densities of sampled Poisson 
#' parameters and weights.
#' 
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with densities for the sampled parameters and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotDens()] for the calling function
".dens.Poisson.Base" <- function(x, dev) {
  K <- x@model@K
  if (.check.grDevice() && dev) {
    dev.new(title = "Densities")
  }
  lambda <- x@par$lambda
  weight <- x@weight
  vars <- cbind(lambda, weight[, seq(1, K - 1)])
  lab.names <- vector("list", 2 * K - 1)
  for (k in seq(1, K)) {
    lab.names[[k]] <- bquote(lambda[.(k)])
  }
  for (k in seq(K + 1, 2 * K - 1)) {
    lab.names[[k]] <- bquote(eta[.(k - K)])
  }
  .symmetric.Dens(vars, lab.names)
}

#' Plot densities of Binomial samples
#' 
#' @description 
#' For internal usage only. This function plots densities of sampled Binomial 
#' parameters and weights.
#' 
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with densities for the sampled parameters and weights.
#' @noRd
#' @importFrom grDevices dev.new
#' @seealso 
#' * [plotDens()] for the calling function
".dens.Binomial.Base" <- function(x, dev) {
  K <- x@model@K
  if (.check.grDevice() && dev) {
    dev.new(title = "Densities")
  }
  p <- x@par$p
  weight <- x@weight
  vars <- cbind(p, weight[, seq(1, K - 1)])
  lab.names <- vector("list", 2 * K - 1)
  for (k in seq(1, K)) {
    lab.names[[k]] <- bquote(p[.(k)])
  }
  for (k in seq(K + 1, 2 * K - 1)) {
    lab.names[[k]] <- bquote(eta[.(k - K)])
  }
  .symmetric.Dens(vars, lab.names)
}

#' Plot densities of exponential samples
#' 
#' @description 
#' For internal usage only. This function plots densities of sampled 
#' exponential parameters and weights.
#' 
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with densities for the sampled parameters and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotDens()] for the calling function
".dens.Exponential.Base" <- function(x, dev) {
  K <- x@model@K
  if (.check.grDevice() && dev) {
    dev.new(title = "Densities")
  }
  lambda <- x@par$lambda
  weight <- x@weight
  vars <- cbind(lambda, weight[, seq(1, K - 1)])
  lab.names <- vector("list", 2 * K - 1)
  for (k in seq(1, K)) {
    lab.names[[k]] <- bquote(lambda[.(k)])
  }
  for (k in seq(K + 1, 2 * K - 1)) {
    lab.names[[k]] <- bquote(eta[.(k - K)])
  }
  .symmetric.Dens(vars, lab.names)
}

#' Plot densities of normal samples
#' 
#' @description 
#' For internal usage only. This function plots densities of sampled normal 
#' parameters and weights.
#' 
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with densities for the sampled parameters and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotDens()] for the calling function
".dens.Normal.Base" <- function(x, dev) {
  K <- x@model@K
  .dens.Normal(x, dev)
  if (K > 1) {
    weight <- x@weight
    weights.lab.names <- vector("list", K)
    for (k in 1:K) {
      weights.lab.names[[k]] <- bquote(eta[.(k)])
    }
    if (K > 1) {
      if (.check.grDevice() && dev) {
        dev.new(title = "Densities Weights")
      }
      .symmetric.Dens(weight, weights.lab.names)
    }
  }
}

#' Plot densities of Student-t samples
#' 
#' @description 
#' For internal usage only. This function plots densities of sampled Student-t 
#' parameters and weights.
#' 
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with densities for the sampled parameters and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotDens()] for the calling function
".dens.Student.Base" <- function(x, dev) {
  K <- x@model@K
  .dens.Student(x, dev)
  if (K > 1) {
    weight <- x@weight
    weight.lab.names <- vector("list", K)
    for (k in 1:K) {
      weight.lab.names[[k]] <- bquote(eta[.(k)])
    }
    if (.check.grDevice() && dev) {
      dev.new(title = "Densities Weights")
    }
    .symmetric.Dens(weight, weight.lab.names)
  }
}

#' Plot densities of multivariate normal samples
#' 
#' @description 
#' For internal usage only. This function plots densities of sampled 
#' multivariate normal parameters and weights.
#' 
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with densities for the sampled parameters and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotDens()] for the calling function
".dens.Normult.Base" <- function(x, dev) {
  K <- x@model@K
  .dens.Normult(x, dev)
  if (K > 1) {
    weight <- x@weight
    weight.lab.names <- vector("list", K)
    for (k in 1:K) {
      weight.lab.names[[k]] <- bquote(eta[.(k)])
    }
    if (.check.grDevice() && dev) {
      dev.new(title = "Densities Weights")
    }
    .symmetric.Dens(weight, weight.lab.names)
  }
}

#' Plot densities of multivariate Student-t samples
#' 
#' @description 
#' For internal usage only. This function plots densities of sampled 
#' multivariate Student-t parameters and weights.
#' 
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with densities for the sampled parameters and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotDens()] for the calling function
".dens.Studmult.Base" <- function(x, dev) {
  K <- x@model@K
  .dens.Studmult(x, dev)
  if (K > 1) {
    weight <- x@weight
    weight.lab.names <- vector("list", K)
    for (k in 1:K) {
      weight.lab.names[[k]] <- bquote(eta[.(k)])
    }
    if (.check.grDevice() && dev) {
      dev.new(title = "Densities Weights")
    }
    .symmetric.Dens(weight, weight.lab.names)
  }
}

### Subseq: Creates a subsequence of an MCMC sample.
#' Generates sub-chains from MCMC samples
#' 
#' @description 
#' For internal usage only. This function generates sub-chains from any 
#' `mcmcoutput` object by defining an `index` array specifying how extraction 
#' of sub-samples should be performed. Has errors for some `mcmcoutput` 
#' sub-classes.
#' 
#' @param obj An `mcmcoutput` object containing all MCMC samples.
#' @param index An array specifying the extraction of sub-samples.
#' @return An `mcmcoutput` object containing sub-chains.
#' @noRd
#' 
#' @seealso 
#' * [subseq()] for the calling method
".subseq.Base" <- function(obj, index) {
  M <- dim(obj@weight)[1]
  K <- dim(obj@weight)[2]
  newM <- sum(index)
  obj@log$cdpost <- array(obj@log$cdpost[index],
    dim = c(newM, 1)
  )
  obj@weight <- obj@weight[index, ]
  obj@entropy <- array(obj@entropy[index],
    dim = c(newM, 1)
  )
  obj@ST <- array(obj@ST[index],
    dim = c(newM, 1)
  )
  ## Check which S stay ##
  storeS <- ifelse(!all(is.na(obj@S)), dim(obj@S)[2], 0)
  if (storeS != 0) {
    ms <- M - storeS
    index.S <- index[(ms + 1):M]
    N <- dim(obj@S)[1]
    if (any(index.S)) {
      obj@S <- array(obj@S[, index.S], dim = c(N, storeS))
    } else {
      obj@S <- as.array(NA)
    }
  }
  obj@NK <- obj@NK[index, ]
  return(obj)
}

#' Swaps elements in MCMC samples.
#' 
#' @description 
#' For internal usage only. This function swaps elements for each row in an 
#' `mcmcoutput` object with MCMC samples.
#' 
#' @param obj An `mcmcoutput` object containing all MCMC samples.
#' @param index An array specifying the element swapping.
#' @return An `mcmcoutput` object with swapped elements.
#' @noRd
#' 
#' @seealso 
#' * [swapElements()] for the calling method
".swapElements.Base" <- function(obj, index) {
  ## Rcpp::export 'swap_cc()'
  obj@weight <- swap_cc(obj@weight, index)
  ## Rcpp::export 'swapInd_cc()'
  M <- obj@M
  K <- ncol(index)
  storeS <- ifelse(!all(is.na(obj@S)), dim(obj@S)[2], 0)
  if (storeS != 0) {
    index.S <- matrix(index[(M - storeS + 1):M, ],
      ncol = K, byrow = TRUE
    )
    obj@S <- swapInd_cc(obj@S, index.S)
  }
  ## Rcpp::export 'swapST_cc()'
  obj@ST <- swapST_cc(obj@ST, index)
  ## Rcpp::export 'swap_cc()'
  obj@NK <- swapInteger_cc(obj@NK, index)
  return(obj)
}
