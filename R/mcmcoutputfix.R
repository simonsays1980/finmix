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



#' Finmix `mcmcoutputfix` class
#' 
#' @description
#' This class defines the basic slots for the MCMC sampling output for a 
#' fixed indicator model. 
#' 
#' @slot M An integer defining the number of iterations in MCMC sampling. 
#' @slot burnin An integer defining the number of iterations in the burn-in 
#'   phase of MCMC sampling. These number of sampling steps are not stored 
#'   in the output.
#' @slot ranperm A logical indicating, if MCMC sampling has been performed 
#'   with random permutations of components.
#' @slot par A named list containing the sampled component parameters. 
#' @slot log A named list containing the values of the mixture log-likelihood, 
#'   mixture prior log-likelihood, and the complete data posterior 
#'   log-likelihood.
#' @slot model The `model` object that specifies the finite mixture model for 
#'   whcih MCMC sampling has been performed. 
#' @slot prior The `prior` object defining the prior distributions for the 
#'   component parameters that has been used in MCMC sampling.
#' @exportClass mcmcoutputfix
#' @rdname mcmcoutputfix-class
#' @keywords internal
.mcmcoutputfix <- setClass("mcmcoutputfix",
  representation(
    M = "integer",
    burnin = "integer",
    ranperm = "logical",
    par = "list",
    log = "list",
    model = "model",
    prior = "prior"
  ),
  validity = function(object) {
    ## else: OK
    TRUE
  },
  prototype(
    M = integer(),
    burnin = integer(),
    ranperm = logical(),
    par = list(),
    log = list(),
    model = model(),
    prior = prior()
  )
)

#' Shows a summary of an `mcmcoutputfix` object.
#' 
#' @description
#' Calling [show()] on an `mcmcoutputfix` object gives an overview 
#' of the `mcmcoutputfix` object.
#' 
#' @param object An `mcmcoutputfix` object.
#' @returns A console output listing the slots and summary information about
#'   each of them. 
#' @exportMethod show
#' @noRd
setMethod(
  "show", "mcmcoutputfix",
  function(object) {
    cat("Object 'mcmcoutputfix'\n")
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
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2, 
#'                  indicfix = TRUE)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc(storepost = FALSE)
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Do not use a hierarchical prior.
#' setHier(f_prior) <- FALSE
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
    x = "mcmcoutputfix",
    dev = "ANY",
    lik = "ANY",
    col = "ANY"
  ),
  function(x, dev = TRUE, lik = 1, col = FALSE, ...) {
    dist <- x@model@dist
    if (lik %in% c(0, 1)) {
      if (dist == "poisson") {
        .traces.Poisson(x, dev)
      } else if (dist == "binomial") {
        .traces.Binomial(x, dev)
      } else if (dist == "exponential") {
        .traces.Exponential(x, dev)
      } else if (dist == "normal") {
        .traces.Normal(x, dev)
      } else if (dist == "student") {
        .traces.Student(x, dev)
      } else if (dist == "normult") {
        .traces.Normult(x, dev, col)
      } else if (dist == "studmult") {
        .traces.Studmult(x, dev, col)
      }
    }
    if (lik %in% c(1, 2)) {
      ## log ##
      .traces.Log(x, dev, col)
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
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2, 
#'                  indicfix = TRUE)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc(storepost = FALSE)
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' setHier(f_prior) <- FALSE
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
    x = "mcmcoutputfix",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    dist <- x@model@dist
    if (dist == "poisson") {
      .hist.Poisson(x, dev)
    } else if (dist == "binomial") {
      .hist.Binomial(x, dev)
    } else if (dist == "exponential") {
      .hist.Exponential(x, dev)
    } else if (dist == "normal") {
      .hist.Normal(x, dev)
    } else if (dist == "student") {
      .hist.Student(x, dev)
    } else if (dist == "normult") {
      .hist.Normult(x, dev)
    } else if (dist == "studmult") {
      .hist.Studmult(x, dev)
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
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2, 
#'                  indicfix = TRUE)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc(storepost = FALSE)
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Do not use a hierarchical prior.
#' setHier(f_prior) <- FALSE
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
    x = "mcmcoutputfix",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    dist <- x@model@dist
    if (dist == "poisson") {
      .dens.Poisson(x, dev)
    } else if (dist == "binomial") {
      .dens.Binomial(x, dev)
    } else if (dist == "exponential") {
      .dens.Exponential(x, dev)
    } else if (dist == "normal") {
      .dens.Normal(x, dev)
    } else if (dist == "student") {
      .dens.Student(x, dev)
    } else if (dist == "normult") {
      .dens.Normult(x, dev)
    } else if (dist == "studmult") {
      .dens.Studmult(x, dev)
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
#' @return Point process of the MCMC samples.
#' @exportMethod plotPointProc
#' @noRd
#' 
#' @examples 
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2, 
#'                  indicfix = TRUE)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc(storepost = FALSE)
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Do not use a hierarchical prior.
#' setHier(f_prior) <- FALSE
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
    x = "mcmcoutputfix",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    dist <- x@model@dist
    if (dist == "poisson") {
      .pointproc.Poisson(x, dev)
    } else if (dist == "binomial") {
      .pointproc.Binomial(x, dev)
    }
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
#' @return Sampling representation of the MCMC samples.
#' @exportMethod plotSampRep
#' @noRd
#' 
#' @examples 
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
    x = "mcmcoutputfix",
    dev = "ANY"
  ),
  function(x, dev, ...) {
    dist <- x@model@dist
    if (dist == "poisson") {
      .samprep.Poisson(x, dev)
    } else if (dist == "binomial") {
      .samprep.Binomial(x, dev)
    }
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
    x = "mcmcoutputfix",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    dist <- x@model@dist
    if (dist == "poisson") {
      .postdens.Poisson(x, dev)
    } else if (dist == "binomial") {
      .postdens.Binomial(x, dev)
    }
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
#' @exportMethod subseq
#' @noRd
setMethod(
  "subseq", signature(
    object = "mcmcoutputfix",
    index = "array"
  ),
  function(object, index) {
    .subseq.valid.Arg(object, index)
    dist <- object@model@dist
    object@M <- sum(index)
    ## log ##
    object <- .subseq.Log.Fix(object, index)
    ## par ##
    if (dist == "poisson") {
      .subseq.Poisson(object, index)
    } else if (dist == "binomial") {
      .subseq.Binomial(object, index)
    } else if (dist == "exponential") {
      .subseq.Poisson(object, index)
    } else if (dist == "normal") {
      .subseq.Normal(object, index)
    } else if (dist == "student") {
      .subseq.Student(object, index)
    } else if (dist == "normult") {
      .subseq.Normult(object, index)
    } else if (dist == "studmult") {
      .subseq.Studmult(object, index)
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
#' @exportMethod swapElements
#' @noRd 
setMethod(
  "swapElements", signature(
    object = "mcmcoutputfix",
    index = "array"
  ),
  function(object, index) { ## Check arguments, TODO: .validObject ##
    .swapElements.valid.Arg(object, index)
    if (object@model@K == 1) {
      return(object)
    } else {
      dist <- object@model@dist
      if (dist == "poisson") {
        .swapElements.Poisson(object, index)
      } else if (dist == "binomial") {
        .swapElements.Binomial(object, index)
      } else if (dist == "exponential") {
        .swapElements.Exponential(object, index)
      } else if (dist == "normal") {
        .swapElements.Normal(object, index)
      } else if (dist == "student") {
        .swapElements.Student(object, index)
      } else if (dist == "normult") {
        .swapElements.Normult(object, index)
      } else if (dist == "studmult") {
        .swapElements.Studmult(object, index)
      }
    }
  }
)

#' Extracts samples from `mcmcoutput` object of a multivariate Normal mixture
#' 
#' @description
#' This function extracts samples from a multivariate Normal mixture output. 
#' 
#' @param object An `mcmcoutput` object from MCMC sampling of a multivariate 
#'   Normal mixture model. 
#' @param index An numeric indicating which dimension of the multivariate 
#'   mixture should be extracted.
#' @return An object class `mcmcextract` containing all samples of an extracted 
#'   dimension.
#' @noRd
#' @exportMethod extract
setMethod(
  "extract", signature(
    object = "mcmcoutputfix",
    index = "numeric"
  ),
  function(object, index) {
    dist <- object@model@dist
    if (dist == "normult") {
      .extract.Normult(object, as.integer(index))
    }
  }
)

# TODO: Check the return values. It appears that this function does not 
# return anything.
#' Computes multivariate Normal sample moments
#' 
#' @description 
#' Calling `moments()` calculates the sample moments for the samples of a 
#' multivariate Normal mixture model.  
#' 
#' @param object An `mcmcoutputfix` object containing all data from MCMC 
#'   sampling.
#' @return The moments on the samples of a multivariate Normal mixture.
#' @exportMethod moments
#' @noRd
setMethod(
  "moments", signature(object = "mcmcoutputfix"),
  function(object) {
    dist <- object@model@dist
    if (dist == "normult") {
      .moments.Normult.Mcmcoutput(object)
    }
  }
)

## Getters ##
#' Getter method of `mcmcoutput` class.
#' 
#' Returns the `M` slot.
#' 
#' @param object An `mcmcoutput` object.
#' @returns The `M` slot of the `object`.
#' @exportMethod getM
#' @noRd
#' 
#' @examples 
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
#' # Get the slot.
#' getM(f_output)
#' 
#' @seealso 
#' * [mcmcoutput-class] for the class definition
#' * [mixturemcmc()] for performing MCMC sampling
setMethod(
  "getM", "mcmcoutputfix",
  function(object) {
    return(object@M)
  }
)

#' Getter method of `mcmcoutput` class.
#' 
#' Returns the `burnin` slot.
#' 
#' @param object An `mcmcoutput` object.
#' @returns The `burnin` slot of the `object`.
#' @exportMethod getBurnin
#' @noRd
#' 
#' @examples 
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
#' # Get the slot.
#' getBurnin(f_output)
#' 
#' @seealso 
#' * [mcmcoutput-class] for the class definition
#' * [mixturemcmc()] for performing MCMC sampling
setMethod(
  "getBurnin", "mcmcoutputfix",
  function(object) {
    return(object@burnin)
  }
)

#' Getter method of `mcmcoutput` class.
#' 
#' Returns the `ranperm` slot.
#' 
#' @param object An `mcmcoutput` object.
#' @returns The `ranperm` slot of the `object`.
#' @exportMethod getRanperm
#' @noRd
#' 
#' @examples 
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
#' # Get the slot.
#' getRanperm(f_output)
#' 
#' @seealso 
#' * [mcmcoutput-class] for the class definition
#' * [mixturemcmc()] for performing MCMC sampling
setMethod(
  "getRanperm", "mcmcoutputfix",
  function(object) {
    return(object@ranperm)
  }
)

#' Getter method of `mcmcoutput` class.
#' 
#' Returns the `par` slot.
#' 
#' @param object An `mcmcoutput` object.
#' @returns The `par` slot of the `object`.
#' @exportMethod getPar
#' @noRd
#' 
#' @examples 
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
#' # Get the slot.
#' getPar(f_output)
#' 
#' @seealso 
#' * [mcmcoutput-class] for the class definition
#' * [mixturemcmc()] for performing MCMC sampling
setMethod(
  "getPar", "mcmcoutputfix",
  function(object) {
    return(object@par)
  }
)

#' Getter method of `mcmcoutput` class.
#' 
#' Returns the `log` slot.
#' 
#' @param object An `mcmcoutput` object.
#' @returns The `log` slot of the `object`.
#' @exportMethod getLog
#' @noRd
#' 
#' @examples 
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
#' # Get the slot.
#' getLog(f_output)
#' 
#' @seealso 
#' * [mcmcoutput-class] for the class definition
#' * [mixturemcmc()] for performing MCMC sampling
setMethod(
  "getLog", "mcmcoutputfix",
  function(object) {
    return(object@log)
  }
)

#' Getter method of `mcmcoutput` class.
#' 
#' Returns the `model` slot.
#' 
#' @param object An `mcmcoutput` object.
#' @returns The `model` slot of the `object`.
#' @exportMethod getModel
#' @noRd
#' 
#' @examples 
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc(storepost = FALSE)
#' # Complete object slots for consistency. 
#' (f_data ~ f_model ~ f_mcmc) %=% mcmcstart(f_data, f_model, f_mcmc)
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' # Get the slot.
#' getModel(f_output)
#' 
#' @seealso 
#' * [mcmcoutput-class] for the class definition
#' * [mixturemcmc()] for performing MCMC sampling
setMethod(
  "getModel", "mcmcoutputfix",
  function(object) {
    return(object@model)
  }
)

#' Getter method of `mcmcoutput` class.
#' 
#' Returns the `prior` slot.
#' 
#' @param object An `mcmcoutput` object.
#' @returns The `prior` slot of the `object`.
#' @exportMethod getPrior
#' @noRd
#' 
#' @examples 
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
#' # Get the slot.
#' getPrior(f_output)
#' 
#' @seealso 
#' * [mcmcoutput-class] for the class definition
#' * [mixturemcmc()] for performing MCMC sampling
setMethod(
  "getPrior", "mcmcoutputfix",
  function(object) {
    return(object@prior)
  }
)

## No setters as users are not intended to manipulate 
## this object 

### Private functions
### These functions are not exported

### Plot

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
".traces.Poisson" <- function(x, dev) {
  K <- x@model@K
  trace.n <- K
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
".traces.Binomial" <- function(x, dev) {
  K <- x@model@K
  trace.n <- K
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
".traces.Exponential" <- function(x, dev) {
  K <- x@model@K
  trace.n <- K
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
  axis(1)
  mtext(side = 1, "Iterations", cex = .7, line = 3)
}

#' Plots traces of normal mixture samples
#' 
#' @description 
#' For internal usage only. This function plots the traces for sampled values 
#' from a normal mixture model.
#' 
#' @param x An `mcmcoutput` object containing all samples.
#' @param dev A logical indicating if the plot should be shown by a grapical 
#'   device.
#' @return A plot of the traces of sampled values.
#' @noRd
#' 
#' @seealso 
#' * [plotTraces()] for the calling function
".traces.Normal" <- function(x, dev) {
  K <- x@model@K
  trace.n <- 2 * K
  if (.check.grDevice() && dev) {
    dev.new(title = "Traceplots")
  }
  par(
    mfrow = c(trace.n, 1), mar = c(1, 0, 0, 0),
    oma = c(4, 5, 4, 4)
  )
  mu <- x@par$mu
  sigma <- x@par$sigma
  for (k in 1:K) {
    plot(mu[, k],
      type = "l", axes = F,
      col = "gray20", xlab = "", ylab = ""
    )
    axis(2, las = 2, cex.axis = .7)
    mtext(
      side = 2, las = 2, bquote(mu[k = .(k)]),
      cex = .6, line = 3
    )
  }
  for (k in 1:K) {
    plot(sigma[, k],
      type = "l", axes = F,
      col = "gray30", xlab = "", ylab = ""
    )
    axis(2, las = 2, cex.axis = .7)
    mtext(
      side = 2, las = 2, bquote(sigma[k = .(k)]),
      cex = .6, line = 3
    )
  }
  axis(1)
  mtext(side = 1, "Iterations", cex = .7, line = 3)
}

#' Plots traces of Student-t mixture samples
#' 
#' @description 
#' For internal usage only. This function plots the traces for sampled values 
#' from a Student-t mixture model.
#' 
#' @param x An `mcmcoutput` object containing all samples.
#' @param dev A logical indicating if the plot should be shown by a grapical 
#'   device.
#' @return A plot of the traces of sampled values.
#' @noRd
#' 
#' @seealso 
#' * [plotTraces()] for the calling function
".traces.Student" <- function(x, dev) {
  K <- x@model@K
  trace.n <- 3 * K
  if (.check.grDevice() && dev) {
    dev.new(title = "Traceplots")
  }
  par(
    mfrow = c(trace.n, 1), mar = c(1, 0, 0, 0),
    oma = c(4, 5, 4, 4)
  )
  mu <- x@par$mu
  sigma <- x@par$sigma
  df <- x@par$df
  for (k in 1:K) {
    plot(mu[, k],
      type = "l", axes = F,
      col = "gray20", xlab = "", ylab = ""
    )
    axis(2, las = 2, cex.axis = .7)
    mtext(
      side = 2, las = 2, bquote(mu[k = .(k)]),
      cex = .6, line = 3
    )
  }
  for (k in 1:K) {
    plot(sigma[, k],
      type = "l", axes = F,
      col = "gray30", xlab = "", ylab = ""
    )
    axis(2, las = 2, cex.axis = .7)
    mtext(
      side = 2, las = 2, bquote(sigma[k = .(k)]),
      cex = .6, line = 3
    )
  }
  for (k in 1:K) {
    plot(df[, k],
      type = "l", axes = F,
      col = "gray40", xlab = "", ylab = ""
    )
    axis(2, las = 2, cex.axis = .7)
    mtext(
      side = 2, las = 2, bquote(nu[k = .(k)]),
      cex = .6, line = 3
    )
  }
  axis(1)
  mtext(side = 1, "Iterations", cex = .7, line = 3)
}

#' Plots traces of multivariate normal mixture samples
#' 
#' @description 
#' For internal usage only. This function plots the traces for sampled values 
#' from a multivariate normal mixture model.
#' 
#' @param x An `mcmcoutput` object containing all samples.
#' @param dev A logical indicating if the plot should be shown by a grapical 
#'   device.
#' @return A plot of the traces of sampled values.
#' @noRd
#' @importFrom grDevices rainbow gray.colors
#' @seealso 
#' * [plotTraces()] for the calling function
".traces.Normult" <- function(x, dev, col) {
  K <- x@model@K
  r <- x@model@r
  if (.check.grDevice() && dev) {
    dev.new(title = "Traceplots")
  }
  trace.n <- r + 2
  par(
    mfrow = c(trace.n, 1), mar = c(1, 2, 0, 0),
    oma = c(4, 5, 2, 4)
  )
  mu <- x@par$mu
  sigma <- x@par$sigma
  if (col) {
    cscale <- rainbow(K, start = 0.5, end = 0)
  } else {
    cscale <- gray.colors(K, start = 0.5, end = 0.15)
  }
  for (rr in 1:r) {
    mmax <- max(mu[, rr, ])
    mmin <- min(mu[, rr, ])
    plot(mu[, rr, 1],
      type = "l", axes = F,
      col = cscale[1], xlab = "", ylab = "",
      ylim = c(mmin, mmax + 0.3 * (mmax - mmin))
    )
    for (k in 2:K) {
      lines(mu[, rr, k], col = cscale[k])
    }
    axis(2, las = 2, cex.axis = .7)
    mtext(
      side = 2, las = 2, bquote(mu[rr = .(rr)]),
      cex = .6, line = 3
    )
    if (rr == 1) {
      name <- vector("character", K)
      for (k in 1:K) {
        name[k] <- paste("k = ", k, sep = "")
      }
      legend("top",
        legend = name, col = cscale, horiz = TRUE,
        lty = 1
      )
    }
  }
  sigma.tr <- array(numeric(), dim = c(x@M, K))
  sigma.det <- array(numeric(), dim = c(x@M, K))
  for (k in 1:K) {
    sigma.tr[, k] <- sapply(
      seq(1, x@M),
      function(i) sum(diag(qinmatr(sigma[i, , k])))
    )
    sigma.det[, k] <- sapply(
      seq(1, x@M),
      function(i) log(det(qinmatr(sigma[i, , k])))
    )
  }
  # Sigma traces
  mmax <- max(sigma.tr)
  mmin <- min(sigma.tr)
  plot(sigma.tr[, 1],
    type = "l", axes = F,
    col = cscale[1], xlab = "", ylab = "",
    ylim = c(mmin, mmax)
  )
  for (k in 2:K) {
    lines(sigma.tr[, k], col = cscale[k])
  }
  axis(2, las = 2, cex.axis = .7)
  mtext(
    side = 2, las = 2, bquote(tr(Sigma)),
    cex = .6, line = 3
  )

  # Sigma logdets
  mmax <- max(sigma.det)
  mmin <- min(sigma.det)
  plot(sigma.det[, 1],
    type = "l", axes = F,
    col = cscale[1], xlab = "", ylab = "",
    ylim = c(mmin, mmax)
  )
  for (k in 2:K) {
    lines(sigma.det[, k], col = cscale[k])
  }
  axis(2, las = 2, cex.axis = .7)
  mtext(
    side = 2, las = 2, bquote(log(det(Sigma))),
    cex = .6, line = 3
  )

  axis(1)
  mtext(side = 1, "Iterations", cex = .7, line = 3)

  # Get moments
  moms <- moments_cc(x)
  for (rr in 1:r) {
    if (.check.grDevice() && dev) {
      dev.new(title = paste("Traceplots Feature ", rr, sep = ""))
    }
    par(
      mfrow = c(2, 2), mar = c(4, 4, 0.5, 0.5),
      oma = c(1.5, 2, 1, 1)
    )
    # Mu
    plot(moms$mean[, rr],
      type = "l", axes = F,
      xlab = "", ylab = "", col = cscale[K]
    )
    axis(2, las = 2, cex.axis = .7)
    mtext(
      side = 2, las = 2, bquote(mu), cex = .6,
      line = 3
    )
    axis(1)
    mtext(side = 1, "Iterations", cex = .7, line = 3)
    # Variance
    plot(moms$var[, rr],
      type = "l", axes = F,
      xlab = "", ylab = "", col = cscale[K]
    )
    axis(2, las = 2, cex.axis = .7)
    mtext(
      side = 2, las = 2, bquote(sigma), cex = .6,
      line = 3
    )
    axis(1)
    mtext(side = 1, "Iterations", cex = .7, line = 3)
    # Skewness
    plot(moms$skewness[, rr],
      type = "l", axes = F,
      xlab = "", ylab = "", col = cscale[K]
    )
    axis(2, las = 2, cex.axis = .7)
    mtext(
      side = 2, "Skewness", cex = .6,
      line = 3
    )
    axis(1)
    mtext(side = 1, "Iterations", cex = .7, line = 3)
    # Kurtosis
    plot(moms$kurtosis[, rr],
      type = "l", axes = F,
      xlab = "", ylab = "", col = cscale[K]
    )
    axis(2, las = 2, cex.axis = .7)
    mtext(
      side = 2, "Kurtosis", cex = .6,
      line = 3
    )
    axis(1)
    mtext(side = 1, "Iterations", cex = .7, line = 3)
  }
}

#' Plots traces of multivariate Student-t mixture samples
#' 
#' @description 
#' For internal usage only. This function plots the traces for sampled values 
#' from a multivariate Student-t mixture model.
#' 
#' @param x An `mcmcoutput` object containing all samples.
#' @param dev A logical indicating if the plot should be shown by a grapical 
#'   device.
#' @return A plot of the traces of sampled values.
#' @noRd
#' 
#' @seealso 
#' * [plotTraces()] for the calling function
".traces.Studmult" <- function(x, dev, col) {
  K <- x@model@K
  r <- x@model@r
  if (.check.grDevice() && dev) {
    dev.new(title = "Traceplots")
  }
  trace.n <- r + 2
  par(
    mfrow = c(trace.n, 1), mar = c(1, 2, 0, 0),
    oma = c(4, 5, 4, 4)
  )
  mu <- x@par$mu
  sigma <- x@par$sigma
  if (col) {
    cscale <- rainbow(K, start = 0.5, end = 0)
  } else {
    cscale <- gray.colors(K, start = 0.5, end = 0.15)
  }
  for (rr in 1:r) {
    mmax <- max(mu[, rr, ])
    mmin <- min(mu[, rr, ])
    plot(mu[, rr, 1],
      type = "l", axes = F,
      col = cscale[1], xlab = "", ylab = "",
      ylim = c(mmin, mmax + 0.3 * (mmax - mmin))
    )
    for (k in 2:K) {
      lines(mu[, rr, k], col = cscale[k])
    }
    axis(2, las = 2, cex.axis = .7)
    mtext(
      side = 2, las = 2, bquote(mu[rr = .(rr)]),
      cex = .6, line = 3
    )
    if (rr == 1) {
      name <- vector("character", K)
      for (k in 1:K) {
        name[k] <- paste("k = ", k, sep = "")
      }
      legend("top",
        legend = name, col = cscale, horiz = TRUE,
        lty = 1
      )
    }
  }
  sigma.tr <- array(numeric(), dim = c(x@M, K))
  sigma.det <- array(numeric(), dim = c(x@M, K))
  for (k in 1:K) {
    sigma.tr[, k] <- sapply(
      seq(1, x@M),
      function(i) sum(diag(qinmatr(sigma[i, , k])))
    )
    sigma.det[, k] <- sapply(
      seq(1, x@M),
      function(i) log(det(qinmatr(sigma[i, , k])))
    )
  }
  # Sigma traces
  mmax <- max(sigma.tr)
  mmin <- min(sigma.tr)
  plot(sigma.tr[, 1],
    type = "l", axes = F,
    col = cscale[1], xlab = "", ylab = "",
    ylim = c(mmin, mmax)
  )
  for (k in 2:K) {
    lines(sigma.tr[, k], col = cscale[k])
  }
  axis(2, las = 2, cex.axis = .7)
  mtext(
    side = 2, las = 2, bquote(tr(Sigma)),
    cex = .6, line = 3
  )

  # Sigma logdets
  mmax <- max(sigma.det)
  mmin <- min(sigma.det)
  plot(sigma.det[, 1],
    type = "l", axes = F,
    col = cscale[1], xlab = "", ylab = "",
    ylim = c(mmin, mmax)
  )
  for (k in 2:K) {
    lines(sigma.det[, k], col = cscale[k])
  }
  axis(2, las = 2, cex.axis = .7)
  mtext(
    side = 2, las = 2, bquote(log(det(Sigma))),
    cex = .6, line = 3
  )

  axis(1)
  mtext(side = 1, "Iterations", cex = .7, line = 3)

  # Degrees of freedom
  if (.check.grDevice() && dev) {
    dev.new(title = "Traceplots Degrees of Freedom")
  }
  degf <- x@par$df
  par(
    mfrow = c(K, 1), mar = c(1, 2, 0, 0),
    oma = c(4, 5, 4, 4)
  )
  for (k in 1:K) {
    plot(degf[, k],
      type = "l", axes = F,
      col = cscale[K], xlab = "", ylab = ""
    )
    axis(2, las = 2, cex.axis = .7)
    mtext(
      side = 2, las = 2, bquote(nu[k = .(k)]),
      cex = .6, line = 3
    )
  }
  axis(1)
  mtext(side = 1, "Iterations", cex = .7, line = 3)
  # Get moments
  moms <- moments_cc(x)
  for (rr in 1:r) {
    if (.check.grDevice() && dev) {
      dev.new(title = paste("Traceplots Feature ", rr, sep = ""))
    }
    par(
      mfrow = c(2, 2), mar = c(4, 4, 0.5, 0.5),
      oma = c(1.5, 2, 1, 1)
    )
    # Mu
    plot(moms$mean[, rr],
      type = "l", axes = F,
      xlab = "", ylab = "", col = cscale[K]
    )
    axis(2, las = 2, cex.axis = .7)
    mtext(
      side = 2, las = 2, bquote(mu), cex = .6,
      line = 3
    )
    axis(1)
    mtext(side = 1, "Iterations", cex = .7, line = 3)
    # Variance
    plot(moms$var[, rr],
      type = "l", axes = F,
      xlab = "", ylab = "", col = cscale[K]
    )
    axis(2, las = 2, cex.axis = .7)
    mtext(
      side = 2, las = 2, bquote(sigma), cex = .6,
      line = 3
    )
    axis(1)
    mtext(side = 1, "Iterations", cex = .7, line = 3)
    # Skewness
    plot(moms$skewness[, rr],
      type = "l", axes = F,
      xlab = "", ylab = "", col = cscale[K]
    )
    axis(2, las = 2, cex.axis = .7)
    mtext(
      side = 2, "Skewness", cex = .6,
      line = 3
    )
    axis(1)
    mtext(side = 1, "Iterations", cex = .7, line = 3)
    # Kurtosis
    plot(moms$kurtosis[, rr],
      type = "l", axes = F,
      xlab = "", ylab = "", col = cscale[K]
    )
    axis(2, las = 2, cex.axis = .7)
    mtext(
      side = 2, "Kurtosis", cex = .6,
      line = 3
    )
    axis(1)
    mtext(side = 1, "Iterations", cex = .7, line = 3)
  }
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
".traces.Log" <- function(x, dev, col) {
  if (.check.grDevice() && dev) {
    dev.new(title = "Log Likelihood Traceplots")
  }
  if (col) {
    cscale <- rainbow(3, start = 0.5, end = 0)
  } else {
    cscale <- gray.colors(3, start = 0.5, end = 0.15)
  }
  par(
    mfrow = c(2, 1), mar = c(1, 0, 0, 0),
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
  axis(1)
  mtext(side = 1, "Iterations", cex = 0.7, line = 3)
}

### Plot Histogramms

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
".hist.Poisson" <- function(x, dev) {
  K <- x@model@K
  if (.check.grDevice() && dev) {
    dev.new(title = "Histograms")
  }
  lambda <- x@par$lambda
  if (K == 1) {
    .symmetric.Hist(lambda, list(bquote(lambda)))
  } else {
    lab.names <- vector("list", K)
    for (k in 1:K) {
      lab.names[[k]] <- bquote(lambda[.(k)])
    }
    .symmetric.Hist(lambda, lab.names)
  }
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
".hist.Binomial" <- function(x, dev) {
  K <- x@model@K
  if (.check.grDevice() && dev) {
    dev.new(title = "Histograms")
  }
  p <- x@par$p
  if (K == 1) {
    .symmetric.Hist(p, list(bquote(p)))
  } else {
    lab.names <- vector("list", K)
    for (k in 1:K) {
      lab.names[[k]] <- bquote(p[.(k)])
    }
    .symmetric.Hist(p, lab.names)
  }
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
".hist.Exponential" <- function(x, dev) {
  K <- x@model@K
  if (.check.grDevice() && dev) {
    dev.new(title = "Histograms")
  }
  lambda <- x@par$lambda
  if (K == 1) {
    .symmetric.Hist(lambda, list(bquote(lambda)))
  } else {
    lab.names <- vector("list", K)
    for (k in 1:K) {
      lab.names[[k]] <- bquote(lambda[.(k)])
    }
    .symmetric.Hist(lambda, lab.names)
  }
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
".hist.Normal" <- function(x, dev) {
  K <- x@model@K
  mu <- x@par$mu
  sigma <- x@par$sigma
  if (K == 1) {
    if (.check.grDevice() && dev) {
      dev.new(title = "Histogram Mu")
    }
    .symmetric.Hist(mu, list(bquote(mu)))
    if (.check.grDevice() && dev) {
      dev.new(title = "Histogram Sigma")
    }
    .symmetric.Hist(sigma, list(bquote(sigma)))
  } else {
    mu.lab.names <- vector("list", K)
    sigma.lab.names <- vector("list", K)
    for (k in 1:K) {
      mu.lab.names[[k]] <- bquote(mu[.(k)])
      sigma.lab.names[[k]] <- bquote(sigma[.(k)])
    }
    if (.check.grDevice() && dev) {
      dev.new(title = "Histograms Mu")
    }
    .symmetric.Hist(mu, mu.lab.names)
    if (.check.grDevice() && dev) {
      dev.new(title = "Histograms Sigma")
    }
    .symmetric.Hist(sigma, sigma.lab.names)
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
".hist.Student" <- function(x, dev) {
  K <- x@model@K
  mu <- x@par$mu
  sigma <- x@par$sigma
  degf <- x@par$df
  if (K == 1) {
    if (.check.grDevice() && dev) {
      dev.new(title = "Histogram Mu")
    }
    .symmetric.Hist(mu, list(bquote(mu)))
    if (.check.grDevice() && dev) {
      dev.new(title = "Histogram Sigma")
    }
    .symmetric.Hist(sigma, list(bquote(sigma)))
    if (.check.grDevice() && dev) {
      dev.new(title = "Histogram Degrees of Freedom")
    }
    .symmetric.Hist(degf, list(bquote(nu)))
  } else {
    mu.lab.names <- vector("list", K)
    sigma.lab.names <- vector("list", K)
    degf.lab.names <- vector("list", K)
    for (k in 1:K) {
      mu.lab.names[[k]] <- bquote(mu[.(k)])
      sigma.lab.names[[k]] <- bquote(sigma[.(k)])
      degf.lab.names[[k]] <- bquote(nu[.(k)])
    }
    if (.check.grDevice() && dev) {
      dev.new(title = "Histograms Mu")
    }
    .symmetric.Hist(mu, mu.lab.names)
    if (.check.grDevice() && dev) {
      dev.new(title = "Histograms Sigma")
    }
    .symmetric.Hist(sigma, sigma.lab.names)
    if (.check.grDevice() && dev) {
      dev.new(title = "Histograms Degrees of Freedom")
    }
    .symmetric.Hist(degf, degf.lab.names)
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
".hist.Normult" <- function(x, dev) {
  K <- x@model@K
  r <- x@model@r
  mu <- x@par$mu
  sigma <- x@par$sigma
  for (rr in 1:r) {
    if (K == 1) {
      if (.check.grDevice() & dev) {
        dev.new(title = paste("Histograms Feature ", rr,
          " Mu",
          sep = ""
        ))
      }
      .symmetric.Hist(mu[, rr, ], list(bquote(mu)))
      if (.check.grDevice() & dev) {
        dev.new(title = paste("Histograms Feature ", rr,
          " Sigma",
          sep = ""
        ))
      }
      .symmetric.Hist(sigma[, rr, ], list(bquote(sigma)))
    } else {
      mu.lab.names <- vector("list", K)
      sigma.lab.names <- vector("list", K)
      for (k in 1:K) {
        mu.lab.names[[k]] <- bquote(mu[.(k)])
        sigma.lab.names[[k]] <- bquote(sigma[.(k)])
      }
      if (.check.grDevice() & dev) {
        dev.new(title = paste("Histograms Feature ", rr,
          " Mu",
          sep = ""
        ))
      }
      .symmetric.Hist(mu[, rr, ], mu.lab.names)
      if (.check.grDevice() & dev) {
        dev.new(title = paste("Histograms Feature ", rr,
          " Sigma",
          sep = ""
        ))
      }
      .symmetric.Hist(sigma[, rr, ], sigma.lab.names)
    }
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
".hist.Studmult" <- function(x, dev) {
  K <- x@model@K
  r <- x@model@r
  mu <- x@par$mu
  sigma <- x@par$sigma
  degf <- x@par$df
  for (rr in 1:r) {
    if (K == 1) {
      if (.check.grDevice() & dev) {
        dev.new(title = paste("Histograms Feature ", rr,
          " Mu",
          sep = ""
        ))
      }
      .symmetric.Hist(mu[, rr, ], list(bquote(mu)))
      if (.check.grDevice() & dev) {
        dev.new(title = paste("Histograms Feature ", rr,
          " Sigma",
          sep = ""
        ))
      }
      .symmetric.Hist(sigma[, rr, ], list(bquote(sigma)))
    } else {
      mu.lab.names <- vector("list", K)
      sigma.lab.names <- vector("list", K)
      for (k in 1:K) {
        mu.lab.names[[k]] <- bquote(mu[.(k)])
        sigma.lab.names[[k]] <- bquote(sigma[.(k)])
      }
      if (.check.grDevice() & dev) {
        dev.new(title = paste("Histograms Feature ", rr,
          " Mu",
          sep = ""
        ))
      }
      .symmetric.Hist(mu[, rr, ], mu.lab.names)
      if (.check.grDevice() & dev) {
        dev.new(title = paste("Histograms Feature ", rr,
          " Sigma",
          sep = ""
        ))
      }
      .symmetric.Hist(sigma[, rr, ], sigma.lab.names)
    }
  }
  if (K == 1) {
    if (.check.grDevice() & dev) {
      dev.new(title = paste("Histograms Feature ", rr,
        " Mu",
        sep = ""
      ))
    }
    .symmetric.Hist(degf[, rr, ], list(bquote(nu)))
  } else {
    degf.lab.names <- vector("list", K)
    for (k in 1:K) {
      degf.lab.names[[k]] <- bquote(nu[.(k)])
    }
    if (.check.grDevice() & dev) {
      dev.new(title = paste("Histograms Feature ", rr,
        " Sigma",
        sep = ""
      ))
    }
    .symmetric.Hist(degf[, rr, ], degf.lab.names)
  }
}

### Plot Densities

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
".dens.Poisson" <- function(x, dev) {
  K <- x@model@K
  if (.check.grDevice() && dev) {
    dev.new(title = "Densities")
  }
  lambda <- x@par$lambda
  if (K == 1) {
    .symmetric.Dens(lambda, list(bquote(lambda)))
  } else {
    lab.names <- vector("list", K)
    for (k in seq(1, K)) {
      lab.names[[k]] <- bquote(lambda[.(k)])
    }
    .symmetric.Dens(lambda, lab.names)
  }
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
#' 
#' @seealso 
#' * [plotDens()] for the calling function
".dens.Binomial" <- function(x, dev) {
  K <- x@model@K
  if (.check.grDevice() && dev) {
    dev.new(title = "Densities")
  }
  p <- x@par$p
  if (K == 1) {
    .symmetric.Dens(p, list(bquote(p)))
  } else {
    lab.names <- vector("list", K)
    for (k in seq(1, K)) {
      lab.names[[k]] <- bquote(p[.(k)])
    }
    .symmetric.Dens(p, lab.names)
  }
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
".dens.Exponential" <- function(x, dev) {
  K <- x@model@K
  if (.check.grDevice() && dev) {
    dev.new(title = "Densities")
  }
  lambda <- x@par$lambda
  if (K == 1) {
    .symmetric.Dens(lambda, list(bquote(lambda)))
  } else {
    lab.names <- vector("list", K)
    for (k in 1:K) {
      lab.names[[k]] <- bquote(lambda[.(k)])
    }
    .symmetric.Dens(lambda, lab.names)
  }
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
".dens.Normal" <- function(x, dev) {
  K <- x@model@K
  mu <- x@par$mu
  sigma <- x@par$sigma
  if (K == 1) {
    if (.check.grDevice() && dev) {
      dev.new(title = "Density Mu")
    }
    .symmetric.Dens(mu, list(bquote(mu)))
    if (.check.grDevice() && dev) {
      dev.new(title = "Density Sigma")
    }
    .symmetric.Dens(sigma, list(bquote(sigma)))
  } else {
    mu.lab.names <- vector("list", K)
    sigma.lab.names <- vector("list", K)
    for (k in 1:K) {
      mu.lab.names[[k]] <- bquote(mu[.(k)])
      sigma.lab.names[[k]] <- bquote(sigma[.(k)])
    }
    if (.check.grDevice() && dev) {
      dev.new(title = "Densities Mu")
    }
    .symmetric.Dens(mu, mu.lab.names)
    if (.check.grDevice() && dev) {
      dev.new(title = "Densities Sigma")
    }
    .symmetric.Dens(sigma, sigma.lab.names)
  }
}

#' Plot densities of Student-t samples with hierarchical priors
#' 
#' @description 
#' For internal usage only. This function plots densities of sampled Student-t
#' parameters and weights when a hierarchical prior had been used.
#' 
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with densities for the sampled parameters and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotDens()] for the calling function
".dens.Student" <- function(x, dev) {
  K <- x@model@K
  mu <- x@par$mu
  sigma <- x@par$sigma
  degf <- x@par$df
  if (K == 1) {
    if (.check.grDevice() && dev) {
      dev.new(title = "Density Mu")
    }
    .symmetric.Dens(mu, list(bquote(mu)))
    if (.check.grDevice() && dev) {
      dev.new(title = "Density Sigma")
    }
    .symmetric.Dens(sigma, list(bquote(sigma)))
    if (.check.grDevice() && dev) {
      dev.new(title = "Density Degrees of Freedom")
    }
    .symmetric.Dens(degf, list(bquote(nu)))
  } else {
    mu.lab.names <- vector("list", K)
    sigma.lab.names <- vector("list", K)
    degf.lab.names <- vector("list", K)
    for (k in 1:K) {
      mu.lab.names[[k]] <- bquote(mu[.(k)])
      sigma.lab.names[[k]] <- bquote(sigma[.(k)])
      degf.lab.names[[k]] <- bquote(nu[.(k)])
    }
    if (.check.grDevice() && dev) {
      dev.new(title = "Densities Mu")
    }
    .symmetric.Dens(mu, mu.lab.names)
    if (.check.grDevice() && dev) {
      dev.new(title = "Densities Sigma")
    }
    .symmetric.Dens(sigma, sigma.lab.names)
    if (.check.grDevice() && dev) {
      dev.new(title = "Densities Degrees of Freedom")
    }
    .symmetric.Dens(degf, degf.lab.names)
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
".dens.Normult" <- function(x, dev) {
  K <- x@model@K
  r <- x@model@r
  mu <- x@par$mu
  sigma <- x@par$sigma
  for (rr in 1:r) {
    if (K == 1) {
      if (.check.grDevice() & dev) {
        dev.new(title = paste("Densities Feature ", rr,
          " Mu",
          sep = ""
        ))
      }
      .symmetric.Dens(mu[, rr, ], list(bquote(mu)))
      if (.check.grDevice() & dev) {
        dev.new(title = paste("Densities Feature ", rr,
          " Sigma",
          sep = ""
        ))
      }
      .symmetric.Dens(sigma[, rr, ], list(bquote(sigma)))
    } else {
      mu.lab.names <- vector("list", K)
      sigma.lab.names <- vector("list", K)
      for (k in 1:K) {
        mu.lab.names[[k]] <- bquote(mu[.(k)])
        sigma.lab.names[[k]] <- bquote(sigma[.(k)])
      }
      if (.check.grDevice() & dev) {
        dev.new(title = paste("Densities Feature ", rr,
          " Mu",
          sep = ""
        ))
      }
      .symmetric.Dens(mu[, rr, ], mu.lab.names)
      if (.check.grDevice() & dev) {
        dev.new(title = paste("Densities Feature ", rr,
          " Sigma",
          sep = ""
        ))
      }
      .symmetric.Dens(sigma[, rr, ], sigma.lab.names)
    }
  }
  if (.check.grDevice() && dev) {
    dev.new(title = "Densities Hyperparameter C")
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
".dens.Studmult" <- function(x, dev) {
  K <- x@model@K
  r <- x@model@r
  mu <- x@par$mu
  sigma <- x@par$sigma
  degf <- x@par$df
  for (rr in 1:r) {
    if (K == 1) {
      if (.check.grDevice() & dev) {
        dev.new(title = paste("Densities Feature ", rr,
          " Mu",
          sep = ""
        ))
      }
      .symmetric.Dens(mu[, rr, ], list(bquote(mu)))
      if (.check.grDevice() & dev) {
        dev.new(title = paste("Densities Feature ", rr,
          " Sigma",
          sep = ""
        ))
      }
      .symmetric.Dens(sigma[, rr, ], list(bquote(sigma)))
    } else {
      mu.lab.names <- vector("list", K)
      sigma.lab.names <- vector("list", K)
      for (k in 1:K) {
        mu.lab.names[[k]] <- bquote(mu[.(k)])
        sigma.lab.names[[k]] <- bquote(sigma[.(k)])
      }
      if (.check.grDevice() & dev) {
        dev.new(title = paste("Densities Feature ", rr,
          " Mu",
          sep = ""
        ))
      }
      .symmetric.Dens(mu[, rr, ], mu.lab.names)
      if (.check.grDevice() & dev) {
        dev.new(title = paste("Densities Feature ", rr,
          " Sigma",
          sep = ""
        ))
      }
      .symmetric.Dens(sigma[, rr, ], sigma.lab.names)
    }
  }
  if (K == 1) {
    if (.check.grDevice() & dev) {
      dev.new(title = "Density Degrees of Freedom")
    }
    .symmetric.Dens(degf[, rr, ], list(bquote(nu)))
  } else {
    degf.lab.names <- vector("list", K)
    for (k in 1:K) {
      degf.lab.names[[k]] <- bquote(nu[.(k)])
    }
    if (.check.grDevice() & dev) {
      dev.new(title = "Densities Degrees of Freedom")
    }
    .symmetric.Dens(degf[, rr, ], degf.lab.names)
  }
}

### Plot Point Processes

#' Plot point processes of Poisson samples
#' 
#' @description 
#' For internal usage only. This function plots the point process of sampled 
#' Poisson parameters and weights against a random normal sample.
#' 
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with the point process for the sampled parameters and weights.
#' @noRd
#' @importFrom stats median
#' @seealso 
#' * [plotPointProc()] for the calling function
".pointproc.Poisson" <- function(x, dev) {
  K <- x@model@K
  M <- x@M
  if (.check.grDevice() && dev) {
    dev.new(title = "Point Process Representation (MCMC)")
  }
  y.grid <- replicate(K, rnorm(M))
  if (median(x@par$lambda) < 1) {
    lambda <- log(x@par$lambda)
  } else {
    lambda <- x@par$lambda
  }
  col.grid <- gray.colors(K,
    start = 0,
    end = 0.5
  )
  legend.names <- vector("list", K)
  for (k in seq(1, K)) {
    legend.names[[k]] <- bquote(lambda[.(k)])
  }
  plot(lambda, y.grid,
    pch = 20, col = col.grid,
    cex = .7, cex.axis = .7, cex.lab = .7,
    main = "", ylab = "", xlab = ""
  )
  mtext(
    side = 1, bquote(lambda), cex = .7,
    cex.lab = .7, line = 3
  )
  legend("topright",
    legend = do.call(
      expression,
      legend.names
    ),
    col = col.grid, fill = col.grid
  )
}

#' Plot point processes of Binomial samples
#' 
#' @description 
#' For internal usage only. This function plots the point process of sampled 
#' Binomial parameters and weights against a random normal sample.
#' 
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with the point process for the sampled parameters and weights.
#' @noRd
#' @importFrom stats rnorm
#' @seealso 
#' * [plotPointProc()] for the calling method
".pointproc.Binomial" <- function(x, dev) {
  K <- x@model@K
  M <- x@M
  if (.check.grDevice() && dev) {
    dev.new(title = "Point Process Representation (MCMC)")
  }
  y.grid <- replicate(K, rnorm(M))
  p <- x@par$p
  col.grid <- gray.colors(K,
    start = 0,
    end = 0.5
  )
  legend.names <- vector("list", K)
  for (k in seq(1, K)) {
    legend.names[[k]] <- bquote(p[.(k)])
  }
  plot(p, y.grid,
    pch = 20, col = col.grid,
    cex = .7, cex.axis = .7, cex.lab = .7,
    main = "", ylab = "", xlab = ""
  )
  mtext(
    side = 1, bquote(p), cex = .7,
    cex.lab = .7, line = 3
  )
  legend("topright",
    legend = do.call(
      expression,
      legend.names
    ),
    col = col.grid, fill = col.grid
  )
}

### Plot sampling representation

#' Plot sampling representation of Poisson samples
#' 
#' @description 
#' For internal usage only. This function plots the sampling representation of 
#' sampled Poisson parameters and weights. Each parameter sample is plotted 
#' against its permuted counterpart.
#' 
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with the sampling representation for the sampled parameters 
#'   and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotSampRep()] for the calling function
".samprep.Poisson" <- function(x, dev) {
  K <- x@model@K
  if (K == 1) {
    warning(paste("Sampling representation is only ",
      "available for mixture models with ",
      "K > 1.",
      sep = ""
    ))
    return(FALSE)
  }
  M <- x@M
  n <- min(2000, x@M)
  n.perm <- choose(K, 2) * factorial(2)
  lambda <- x@par$lambda
  if (.check.grDevice() && dev) {
    dev.new(title = "Sampling Representation (MCMC)")
  }
  comb <- as.matrix(expand.grid(seq(1, K), seq(1, K)))
  comb <- comb[which(comb[, 1] != comb[, 2]), ]
  lambda <- lambda[seq(1, n), ]
  lambda <- matrix(lambda[, comb], nrow = n * n.perm, ncol = 2)
  plot(lambda,
    col = "gray47", cex.lab = .7, cex.axis = .7,
    cex = .7, pch = 20, main = "", xlab = "", ylab = ""
  )
  abline(0, 1, lty = 1)
  mtext(
    side = 1, bquote(lambda), cex = .7, cex.lab = .7,
    line = 3
  )
  mtext(
    side = 2, bquote(lambda), cex = .7, cex.lab = .7,
    line = 3
  )
}

#' Plot sampling representation of Binomial samples
#' 
#' @description 
#' For internal usage only. This function plots the sampling representation of 
#' sampled Binomial parameters and weights. Each parameter sample is plotted 
#' against its permuted counterpart.
#' 
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with the sampling representation for the sampled parameters 
#'   and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotSampRep()] for the calling function
".samprep.Binomial" <- function(x, dev) {
  K <- x@model@K
  if (K == 1) {
    warning(paste("Sampling representation is only ",
      "available for mixture models with ",
      "K > 1.",
      sep = ""
    ))
    return(FALSE)
  }
  M <- x@M
  n <- min(2000, x@M)
  n.perm <- choose(K, 2) * factorial(2)
  p <- x@par$p
  if (.check.grDevice() && dev) {
    dev.new(title = "Sampling Representation")
  }
  comb <- as.matrix(expand.grid(seq(1, K), seq(1, K)))
  comb <- comb[which(comb[, 1] != comb[, 2]), ]
  p <- p[seq(1, n), ]
  p <- matrix(p[, comb], nrow = n * n.perm, ncol = 2)
  plot(p,
    col = "gray47", cex.lab = .7, cex.axis = .7,
    cex = .7, pch = 20, main = "", xlab = "", ylab = ""
  )
  abline(0, 1, lty = 1)
  mtext(
    side = 1, bquote(p), cex = .7, cex.lab = .7,
    line = 3
  )
  mtext(
    side = 2, bquote(p), cex = .7, cex.lab = .7,
    line = 3
  )
}

### Posterior Density

#' Plot posterior density of Poisson samples
#' 
#' @description 
#' For internal usage only. This function plots the posterior density of 
#' sampled Poisson parameters and weights. 
#' 
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with the posterior density for the sampled parameters 
#'   and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotPostdens()] for the calling function
".postdens.Poisson" <- function(x, dev) {
  K <- x@model@K
  if (K != 2) {
    warning(paste("A plot of the posterior density is ",
      "available only for K = 2.",
      sep = ""
    ))
  } else {
    M <- x@M
    n <- min(2000, M)
    lambda <- x@par$lambda
    lambda <- lambda[seq(1, n), ]
    dens <- bkde2D(lambda, bandwidth = c(
      sd(lambda[, 1]),
      sd(lambda[, 2])
    ))
    if (.check.grDevice() && dev) {
      dev.new(title = "Posterior Density Contour Plot (MCMC)")
    }
    contour(dens$x1, dens$x2, dens$fhat,
      cex = .7,
      cex.lab = .7, cex.axis = .7, col = "gray47",
      main = "", xlab = "", ylab = ""
    )
    mtext(
      side = 1, bquote(lambda[1]), cex = .7,
      cex.lab = .7, line = 3
    )
    mtext(
      side = 2, bquote(lambda[2]), cex = .7,
      cex.lab = .7, line = 3
    )
    if (.check.grDevice() && dev) {
      dev.new(title = "Posterior Density Perspective Plot (MCMC)")
    }
    persp(dens$x1, dens$x2, dens$fhat,
      col = "gray65",
      border = "gray47", theta = 55, phi = 30,
      expand = .5, lphi = 180, ltheta = 90,
      r = 40, d = .1, ticktype = "detailed", zlab =
        "Density", xlab = "k = 1", ylab = "k = 2"
    )
  }
}

#' Plot posterior density of Binomial samples
#' 
#' @description 
#' For internal usage only. This function plots the posterior density of 
#' sampled Binomial parameters and weights. 
#' 
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with the posterior density for the sampled parameters 
#'   and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotPostDens()] for the calling function
".postdens.Binomial" <- function(x, dev) {
  K <- x@model@K
  if (K != 2) {
    warning(paste("A plot of the posterior density is ",
      "available only for K = 2.",
      sep = ""
    ))
  } else {
    M <- x@M
    n <- min(2000, M)
    p <- x@par$p
    p <- p[seq(1, n), ]
    dens <- bkde2D(p, bandwidth = c(
      sd(p[, 1]),
      sd(p[, 2])
    ))
    if (.check.grDevice() && dev) {
      dev.new(title = "Posterior Density Contour Plot (MCMC)")
    }
    contour(dens$x1, dens$x2, dens$fhat,
      cex = .7,
      cex.lab = .7, cex.axis = .7, col = "gray47",
      main = "", xlab = "", ylab = ""
    )
    mtext(
      side = 1, bquote(p[1]), cex = .7,
      cex.lab = .7, line = 3
    )
    mtext(
      side = 2, bquote(p[2]), cex = .7,
      cex.lab = .7, line = 3
    )
    if (.check.grDevice() && dev) {
      dev.new(title = "Posterior Density Perspective Plot (MCMC)")
    }
    persp(dens$x1, dens$x2, dens$fhat,
      col = "gray65",
      border = "gray47", theta = 55, phi = 30,
      expand = .5, lphi = 180, ltheta = 90,
      r = 40, d = .1, ticktype = "detailed", zlab =
        "Density", xlab = "k = 1", ylab = "k = 2"
    )
  }
}

### Logic
### Logic subseq: This function is used for each
### distribution type in 'model'. It crreates a subsequence
### for the log-likelihoods.
#' Generates sub-chains from MCMC log-likelihood samples
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
".subseq.Log.Fix" <- function(obj, index) {
  obj@log$mixlik <- matrix(obj@log$mixlik[index,],
    nrow = obj@M, ncol = 1
  )
  obj@log$mixprior <- matrix(obj@log$mixprior[index,],
    nrow = obj@M, ncol = 1
  )
  return(obj)
}

#' Generates sub-chains from Poisson MCMC samples
#' 
#' @description 
#' For internal usage only. This function generates sub-chains from an 
#' `mcmcoutput` object with a Poisson `model` by defining an `index` array 
#' specifying how extraction of sub-samples should be performed. Has errors for 
#' some `mcmcoutput` sub-classes.
#' 
#' @param obj An `mcmcoutput` object containing all MCMC samples.
#' @param index An array specifying the extraction of sub-samples.
#' @return An `mcmcoutput` object containing sub-chains.
#' @noRd
#' 
#' @seealso 
#' * [subseq()] for the calling method
".subseq.Poisson" <- function(obj, index) {
  if (obj@model@K == 1) {
    obj@par$lambda <- matrix(obj@par$lambda[index],
      nrow = obj@M, ncol = 1
    )
  } else {
    obj@par$lambda <- obj@par$lambda[index,]
  }
  return(obj)
}

#' Generates sub-chains from Binomial MCMC samples
#' 
#' @description 
#' For internal usage only. This function generates sub-chains from an 
#' `mcmcoutput` object with a Binomial `model` by defining an `index` array 
#' specifying how extraction of sub-samples should be performed. Has errors for 
#' some `mcmcoutput` sub-classes.
#' 
#' @param obj An `mcmcoutput` object containing all MCMC samples.
#' @param index An array specifying the extraction of sub-samples.
#' @return An `mcmcoutput` object containing sub-chains.
#' @noRd
#' 
#' @seealso 
#' * [subseq()] for the calling method
".subseq.Binomial" <- function(obj, index) {
  if (obj@model@K == 1) {
    obj@par$p <- matrix(obj@par$p[index],
      nrow = obj@M,
      ncol = 1
    )
  } else {
    obj@par$p <- obj@par$p[index, ]
  }
  return(obj)
}

#' Generates sub-chains from normal MCMC samples
#' 
#' @description 
#' For internal usage only. This function generates sub-chains from an 
#' `mcmcoutput` object with a normal `model` by defining an `index` array 
#' specifying how extraction of sub-samples should be performed. Has errors for 
#' some `mcmcoutput` sub-classes.
#' 
#' @param obj An `mcmcoutput` object containing all MCMC samples.
#' @param index An array specifying the extraction of sub-samples.
#' @return An `mcmcoutput` object containing sub-chains.
#' @noRd
#' 
#' @seealso 
#' * [subseq()] for the calling method
".subseq.Normal" <- function(obj, index) {
  if (obj@model@K == 1) {
    obj@par$mu <- matrix(obj@par$mu[index],
      nrow = obj@M,
      ncol = 1
    )
    obj@par$sigma <- matrix(obj@par$mu[index],
      nrow = obj@M,
      ncol = 1
    )
  } else {
    obj@par$mu <- obj@par$mu[index, ]
    obj@par$sigma <- obj@par$sigma[index, ]
  }
  return(obj)
}

#' Generates sub-chains from Student-t MCMC samples
#' 
#' @description 
#' For internal usage only. This function generates sub-chains from an 
#' `mcmcoutput` object with a Student-t `model` by defining an `index` array 
#' specifying how extraction of sub-samples should be performed. Has errors for 
#' some `mcmcoutput` sub-classes.
#' 
#' @param obj An `mcmcoutput` object containing all MCMC samples.
#' @param index An array specifying the extraction of sub-samples.
#' @return An `mcmcoutput` object containing sub-chains.
#' @noRd
#' 
#' @seealso 
#' * [subseq()] for the calling method
".subseq.Student" <- function(obj, index) {
  if (obj@model@K == 1) {
    obj@par$mu <- matrix(obj@par$mu[index],
      nrow = obj@M,
      ncol = 1
    )
    obj@par$sigma <- matrix(obj@par$sigma[index],
      nrow = obj@M,
      ncol = 1
    )
    obj@par$df <- matrix(obj@par$df[index],
      nrow = obj@M,
      ncol = 1
    )
  } else {
    obj@par$mu <- obj@par$mu[index, ]
    obj@par$sigma <- obj@par$sigma[index, ]
    obj@par$df <- obj@par$df[index, ]
  }
  return(obj)
}

#' Generates sub-chains from multivariate Normal MCMC samples
#' 
#' @description 
#' For internal usage only. This function generates sub-chains from an 
#' `mcmcoutput` object with a multivariate Normal `model` by defining an 
#' `index` array specifying how extraction of sub-samples should be performed. 
#' Has errors for some `mcmcoutput` sub-classes.
#' 
#' @param obj An `mcmcoutput` object containing all MCMC samples.
#' @param index An array specifying the extraction of sub-samples.
#' @return An `mcmcoutput` object containing sub-chains.
#' @noRd
#' 
#' @seealso 
#' * [subseq()] for the calling method
".subseq.Normult" <- function(obj, index) {
  if (obj@model@K == 1) {
    obj@par$mu <- matrix(obj@par$mu[index, ],
      nrow = obj@M,
      ncol = 1
    )
    obj@par$sigma <- matrix(obj@par$sigma[index, ],
      nrow = obj@M,
      ncol = 1
    )
    obj@par$sigmainv <- matrix(obj@par$sigmainv[index, ],
      nrow = obj@M,
      ncol = 1
    )
  } else {
    obj@par$mu <- obj@par$mu[index, , ]
    obj@par$sigma <- obj@par$sigma[index, , ]
    obj@par$sigmainv <- obj@par$sigmainv[index, , ]
  }
  return(obj)
}

#' Generates sub-chains from Student-t MCMC samples
#' 
#' @description 
#' For internal usage only. This function generates sub-chains from an 
#' `mcmcoutput` object with a Student-t `model` by defining an `index` array 
#' specifying how extraction of sub-samples should be performed. Has errors for 
#' some `mcmcoutput` sub-classes.
#' 
#' @param obj An `mcmcoutput` object containing all MCMC samples.
#' @param index An array specifying the extraction of sub-samples.
#' @return An `mcmcoutput` object containing sub-chains.
#' @noRd
#' 
#' @seealso 
#' * [subseq()] for the calling method
".subseq.Studmult" <- function(obj, index) {
  if (obj@model@K == 1) {
    obj@par$mu <- obj@par$mu[index, ]
    obj@par$sigma <- obj@par$sigma[index, ]
    obj@par$sigmainv <- obj@par$sigmainv[index, ]
    obj@par$df <- obj@par$df[index]
  } else {
    obj@par$mu <- obj@par$mu[index, , ]
    obj@par$sigma <- obj@par$sigma[index, , ]
    obj@par$sigmainv <- obj@par$sigmainv[index, , ]
    obj@par$df <- obj@par$df[index, ]
  }
  return(obj)
}

#' Swaps elements in Poisson MCMC samples.
#' 
#' @description 
#' For internal usage only. This function swaps elements for each row in an 
#' `mcmcoutput` object with Poisson MCMC samples. Calls the C++-function 
#' `swap_cc()`.
#' 
#' @param obj An `mcmcoutput` object containing all MCMC samples.
#' @param index An array specifying the element swapping.
#' @return An `mcmcoutput` object with swapped elements.
#' @noRd
#' 
#' @seealso 
#' * [swapElements()] for the calling method
".swapElements.Poisson" <- function(obj, index) {
  ## Rcpp::export 'swap_cc'
  obj@par$lambda <- swap_cc(obj@par$lambda, index)
  return(obj)
}

#' Swaps elements in Binomial MCMC samples.
#' 
#' @description 
#' For internal usage only. This function swaps elements for each row in an 
#' `mcmcoutput` object with Binomial MCMC samples. Calls the C++-function 
#' `swap_cc()`.
#' 
#' @param obj An `mcmcoutput` object containing all MCMC samples.
#' @param index An array specifying the element swapping.
#' @return An `mcmcoutput` object with swapped elements.
#' @noRd
#' 
#' @seealso 
#' * [swapElements()] for the calling method
".swapElements.Binomial" <- function(obj, index) {
  ## Rcpp::export 'swap_cc'
  obj@par$p <- swap_cc(obj@par$p, index)
  return(obj)
}

#' Swaps elements in Exponential MCMC samples.
#' 
#' @description 
#' For internal usage only. This function swaps elements for each row in an 
#' `mcmcoutput` object with ExponentialMCMC samples. Calls the C++-function 
#' `swap_cc()`.
#' 
#' @param obj An `mcmcoutput` object containing all MCMC samples.
#' @param index An array specifying the element swapping.
#' @return An `mcmcoutput` object with swapped elements.
#' @noRd
#' 
#' @seealso 
#' * [swapElements()] for the calling method
".swapElements.Exponential" <- function(obj, index) {
  ## Rcpp::export 'swap_cc'
  obj@par$lambda <- swap_cc(obj@par$lambda, index)
  return(obj)
}

#' Swaps elements in Normal MCMC samples.
#' 
#' @description 
#' For internal usage only. This function swaps elements for each row in an 
#' `mcmcoutput` object with Normal MCMC samples. Calls the C++-function 
#' `swap_cc()`.
#' 
#' @param obj An `mcmcoutput` object containing all MCMC samples.
#' @param index An array specifying the element swapping.
#' @return An `mcmcoutput` object with swapped elements.
#' @noRd
#' 
#' @seealso 
#' * [swapElements()] for the calling method
".swapElements.Normal" <- function(obj, index) {
  ## Rcpp::export 'swap_cc'
  obj@par$mu <- swap_cc(obj@par$mu, index)
  obj@par$sigma <- swap_cc(obj@par$sigma, index)
  return(obj)
}

#' Swaps elements in Student-t MCMC samples.
#' 
#' @description 
#' For internal usage only. This function swaps elements for each row in an 
#' `mcmcoutput` object with Student-t MCMC samples. Calls the C++-function 
#' `swap_cc()`.
#' 
#' @param obj An `mcmcoutput` object containing all MCMC samples.
#' @param index An array specifying the element swapping.
#' @return An `mcmcoutput` object with swapped elements.
#' @noRd
#' 
#' @seealso 
#' * [swapElements()] for the calling method
".swapElements.Student" <- function(obj, index) {
  ## Rcpp::export 'swap_cc'
  obj@par$mu <- swap_cc(obj@par$mu, index)
  obj@par$sigma <- swap_cc(obj@par$sigma, index)
  obj@par$df <- swap_cc(obj@par$df, index)
  return(obj)
}

#' Swaps elements in multivariate Normal MCMC samples.
#' 
#' @description 
#' For internal usage only. This function swaps elements for each row in an 
#' `mcmcoutput` object with multivariate Normal MCMC samples. Calls the 
#' C++-function `swap_cc()`.
#' 
#' @param obj An `mcmcoutput` object containing all MCMC samples.
#' @param index An array specifying the element swapping.
#' @return An `mcmcoutput` object with swapped elements.
#' @noRd
#' 
#' @seealso 
#' * [swapElements()] for the calling method
".swapElements.Normult" <- function(obj, index) {
  ## Rcpp::export 'swap_3d_cc'
  obj@par$mu <- swap_3d_cc(obj@par$mu, index)
  obj@par$sigma <- swap_3d_cc(obj@par$sigma, index)
  obj@par$sigmainv <- swap_3d_cc(obj@par$sigma, index)
  return(obj)
}

#' Swaps elements in multivariate Student-t MCMC samples.
#' 
#' @description 
#' For internal usage only. This function swaps elements for each row in an 
#' `mcmcoutput` object with multivariate Student-t MCMC samples. Calls the C++-function 
#' `swap_cc()`.
#' 
#' @param obj An `mcmcoutput` object containing all MCMC samples.
#' @param index An array specifying the element swapping.
#' @return An `mcmcoutput` object with swapped elements.
#' @noRd
#' 
#' @seealso 
#' * [swapElements()] for the calling method
".swapElements.Studmult" <- function(obj, index) {
  ## Rcpp::export 'swap_3d_cc'
  obj@par$mu <- swap_3d_cc(obj@par$mu, index)
  obj@par$sigma <- swap_3d_cc(obj@par$sigma, index)
  obj@par$sigmainv <- swap_3d_cc(obj@par$sigma, index)
  obj@par$df <- swap_cc(obj@par$df, index)
  return(obj)
}

### Validity

#' Checks arguments to `subseq()`
#' 
#' @description 
#' For internal usage only. This functions checks the input arguments to the 
#' [subseq()] method for validity. More specifically, the `index` argument is 
#' checked for the right dimension and type. `index` must be an array of 
#' dimension `M x 1` of logicals indicating which samples should be 
#' extracted.
#' 
#' @param obj An `mcmcoutput` object containing the MCMC sampling that should 
#'   be sub-chained. 
#' @param index An array og logicals indicating which indices should be 
#'   extracted.
#' @return None. If any check does not pass an error is thrown to explain the 
#'   user what to change.
#' @noRd
#' @seealso 
#' * [subseq()] for the calling method
".subseq.valid.Arg" <- function(obj, index) {
  if (dim(index)[1] != obj@M) {
    stop("Argument 'index' has wrong dimension.")
  }
  if (typeof(index) != "logical") {
    stop("Argument 'index' must be of type 'logical'.")
  }
}


#' Checks arguments to `subseq()`
#' 
#' @description 
#' For internal usage only. This functions checks the input arguments to the 
#' [swapElements()] method for validity. More specifically, the `index` 
#' argument is checked for the right dimension and type. `index` must be an 
#' array of dimension `M x K` of integers in the range `1,...,K` indicating 
#' how components should be swapped in each row.
#' 
#' @param obj An `mcmcoutput` object containing the MCMC sampling that should 
#'   be sub-chained. 
#' @param index An array of indicators indicating how components should be 
#'   swapped
#' @return None. If any check does not pass an error is thrown to explain the 
#'   user what to change.
#' @noRd
#' @seealso 
".swapElements.valid.Arg" <- function(obj, index) {
  M <- ifelse(inherits(obj, "mcmcoutputperm"), obj@Mperm, obj@M)
  if (dim(index)[1] != M || dim(index)[2] != obj@model@K) {
    stop("Argument 'index' has wrong dimension.")
  }
  if (typeof(index) != "integer") {
    stop("Argument 'index' must be of type 'integer'.")
  }
  if (!all(index > 0) || any(index > obj@model@K)) {
    stop(paste("Elements in argument 'index' must be greater 0",
      "and must not exceed its number of columns.",
      sep = ""
    ))
  }
}

#' Extract samples from a multivariate Normal `mcmcoutput` object
#' 
#' @description 
#' For internal usage only. This function extracts samples from an `mcmcoutput` 
#' object of a multivariate Normal mixture. 
#' 
#' @param obj An `mcmcoutput` object containing the MCMC samples.
#' @param index An array of logicals indicating which samples should be 
#'   extracted from the MCMC sampling output.
#' @return An `mcmcoutput` object containin the extracted MCMC samples.
#' @noRd
".extract.Normult" <- function(obj, index) {
  dist <- obj@model@dist
  r <- obj@model@r
  K <- obj@model@K
  pars <- sapply(obj@par, function(x) x[index, , ])
  weight <- as.array(obj@model@weight)
  .mcmcextract(
    dist = dist, K = K, r = r, par = pars,
    weight = weight
  )
}

#' Moments for each sample of a multivariate Normal mixture
#' 
#' @description 
#' For internal usage only. This function calculates the moments for extracted 
#' samples from a multivariate Normal mixture model.
#' 
#' @param obj An `mcmcoutput` object containing the MCMC samples.
#' @return An `mcmcextract` object containing the moments.
#' @noRd
#' @seealso 
#' * [mcmcextract][mcmcextract_class]
".moments.Normult.Mcmcoutput" <- function(obj) {
  moments <- array(numeric(), dim = c(obj@M, obj@fdata@r, obj@fdata@r))
  moments <- apply(
    seq(1, obj@M), 1,
    function(i) {
      mm <- extract(obj, i)
      moms <- moments(mm)
    }
  )
}