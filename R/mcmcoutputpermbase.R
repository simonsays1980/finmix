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

#' Finmix `mcmcoutputpermbase` class 
#' 
#' @description 
#' This class defines objects to store the outputs from permuting the MCMC 
#' samples. Due to label switching the sampled component parameters are usually 
#' not assigned to the same component in each iteration. To overcome this issue 
#' the samples are permuted by using a relabeling algorithm (usually K-means) 
#' to reassign parameters. Note that due to assignment of parameters from the 
#' same iteration to the same component, the sample size could shrink. 
#' 
#' This class stores the permuted parameters together with the new MCMC sample 
#' size and the mixture log-likelihood, the prior log-likelihood, and the 
#' complete data posterior log-likelihood.
#' 
#' Note that this class inherits all of its slots from the parent classes.
#' 
#' @exportClass mcmcoutputpermbase
#' @rdname mcmcoutputpermbase-class
#' @seealso 
#' * [mcmcoutputbase-class] for the parent class
#' * [mcmcpermind-class] for the parent class
#' * [mcmcpermute()] for performing permutation of MCMC samples
.mcmcoutputpermbase <- setClass("mcmcoutputpermbase",
  contains = c(
    "mcmcpermind",
    "mcmcoutputbase"
  ),
  validity = function(object) {
    ## else: OK
    TRUE
  }
)

#' Initializer of the `mcmcoutputpermbase` class
#' 
#' @description
#' Only used implicitly. The initializer stores the data into the slots of the 
#' passed-in object.
#' 
#' @param .Object An object: see the "initialize Methods" section in 
#'   [initialize].
#' @param mcmcoutput An `mcmcoutputpermbase` class containing the results from MCMC 
#'   sampling.
#' @param Mperm An integer defining the number of permuted MCMC samples.
#' @param parperm A named list containing the permuted component parameter 
#'   samples from MCMC sampling
#' @param relabel A character specifying the relabeling algorithm used for 
#'   permuting the MCMC samples.
#' @param weightperm An array of dimension `Mperm x K` containing the 
#'   relabeled weight parameters. 
#' @param logperm A named list containing the mixture log-likelihood, the 
#'   prior log-likelihood, and the complete data posterior log-likelihood 
#'   for the permuted MCMC samples.
#' @param entropyperm An `array` of dimension `Mperm x 1` containing the 
#'   entropy for each MCMC permuted draw.
#' @param STperm An `array` of dimension `Mperm x 1` containing all permuted 
#'   MCMC states, for the last observation in slot `@@y` of the `fdata` object 
#'   passed in to [mixturemcmc()] where a state is defined for non-Markov 
#'   models as the last indicator of this observation.  
#' @param Sperm An `array` of dimension `N x storeS` containing the last 
#'   `storeS` permuted indicators. `storeS` is defined in the slot `@@storeS` 
#'   of the `mcmc` object passed into [mixturemcmc()].
#' @param NKperm 
#' 
#' @keywords internal
#' 
#' @seealso 
#' * [Classes_Details] for details of class definitions, and 
#' * [setOldClass] for the relation to S3 classes
setMethod(
  "initialize", "mcmcoutputpermbase",
  function(.Object, mcmcoutput, Mperm = integer(),
           parperm = list(), relabel = character(),
           weightperm = array(), logperm = list(),
           entropyperm = array(), STperm = array(),
           Sperm = array(), NKperm = array()) {
    .Object@M <- mcmcoutput@M
    .Object@burnin <- mcmcoutput@burnin
    .Object@ranperm <- mcmcoutput@ranperm
    .Object@par <- mcmcoutput@par
    .Object@weight <- mcmcoutput@weight
    .Object@log <- mcmcoutput@log
    .Object@ST <- mcmcoutput@ST
    .Object@S <- mcmcoutput@S
    .Object@NK <- mcmcoutput@NK
    .Object@clust <- mcmcoutput@clust
    .Object@model <- mcmcoutput@model
    .Object@prior <- mcmcoutput@prior
    .Object@Mperm <- Mperm
    .Object@parperm <- parperm
    .Object@relabel <- relabel
    .Object@weightperm <- weightperm
    .Object@logperm <- logperm
    .Object@entropyperm <- entropyperm
    .Object@STperm <- STperm
    .Object@Sperm <- Sperm
    .Object@NKperm <- NKperm
    .Object
  }
)

#' Shows a summary of an `mcmcoutputpermbase` object.
#' 
#' @description
#' Calling [show()] on an `mcmcoutputpermbase` object gives an overview 
#' of the `mcmcoutputpermbase` object.
#' 
#' @param object An `mcmcoutputpermbase` object.
#' @returns A console output listing the slots and summary information about
#'   each of them. 
#' @exportMethod show
#' @describeIn mcmcoutputpermbase-class Shows a short summary of the object's 
#'   slots
setMethod(
  "show", "mcmcoutputpermbase",
  function(object) {
    cat("Object 'mcmcoutputperm'\n")
    cat("     class       :", class(object), "\n")
    cat("     M           :", object@M, "\n")
    cat("     burnin      :", object@burnin, "\n")
    cat("     ranperm     :", object@ranperm, "\n")
    cat("     relabel     :", object@relabel, "\n")
    cat(
      "     par         : List of",
      length(object@par), "\n"
    )
    cat(
      "     log         : List of",
      length(object@log), "\n"
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
    cat("     Mperm       :", object@Mperm, "\n")
    cat(
      "     parperm     : List of",
      length(object@parperm), "\n"
    )
    cat(
      "     weightperm  :",
      paste(dim(object@weightperm), collapse = "x"), "\n"
    )
    cat(
      "     logperm     : List of",
      length(object@logperm), "\n"
    )
    cat(
      "     entropyperm :",
      paste(dim(object@entropyperm), collapse = "x"), "\n"
    )
    cat(
      "     STperm      :",
      paste(dim(object@STperm), collapse = "x"), "\n"
    )
    if (!all(is.na(object@Sperm))) {
      cat(
        "     Sperm       :",
        paste(dim(object@Sperm), collapse = "x"), "\n"
      )
    }
    cat(
      "     NKperm      :",
      paste(dim(object@NKperm), collapse = "x"), "\n"
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
#' set to `1`.
#' 
#' If `lik` is set to `0` the parameters of the components and the posterior 
#' parameters are plotted together with `K-1` weights.
#' 
#' @param x An `mcmcoutputpermbase` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param lik An integer indicating, if the log-likelihood traces should be 
#'   plotted (default). If set to `0` the traces for the parameters 
#'   and weights are plotted instead. 
#' @param col A logical indicating, if the plot should be colored.
#' @param ... Further arguments to be passed to the plotting function.
#' @return A plot of the traces of the MCMC samples.
#' @exportMethod plotTraces
#' @describeIn mcmcoutputpermbase-class
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
#' # Do not use a hierarchical prior.
#' setHier(f_prior) <- FALSE
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' f_outputperm <- mcmcpermute(f_output)
#' plotTraces(f_outputperm, lik = 0)
#' }
#' 
#' @seealso 
#' * [mixturemcmc()] for performing MCMC sampling
#' * [mcmcpermute()] for permuting MCMC samples 
#' * [plotHist()] for plotting histograms of sampled values
#' * [plotDens()] for plotting densities of sampled values
#' * [plotSampRep()] for plotting sampling representations of sampled values
#' * [plotPointProc()] for plotting point processes for sampled values
#' * [plotPostDens()] for plotting the posterior density of component parameters
setMethod(
  "plotTraces", signature(
    x = "mcmcoutputpermbase",
    dev = "ANY",
    lik = "ANY",
    col = "ANY"
  ),
  function(x, dev = TRUE, lik = 1, col = FALSE, ...) {
    dist <- x@model@dist
    if (lik %in% c(0, 1)) {
      if (dist == "poisson") {
        .permtraces.Poisson.Base(x, dev)
      } else if (dist == "binomial") {
        .permtraces.Binomial.Base(x, dev)
      } else if (dist == "exponential") {
        .permtraces.Exponential.Base(x, dev)
      } else if (dist == "normal") {
        .permtraces.Normal(x, dev)
        .permtraces.Weights.Base(x, dev, col)
      } else if (dist == "student") {
        .permtraces.Student(x, dev)
        .permtraces.Weights.Base(x, dev, col)
      } else if (dist == "normult") {
        .permtraces.Normult(x, dev, col)
        .permtraces.Weights.Base(x, dev, col)
      } else if (dist == "studmult") {
        .permtraces.Studmult(x, dev, col)
        .permtraces.Weights.Base(x, dev, col)
      }
    }
    if (lik %in% c(1, 2)) {
      ## log ##
      .permtraces.Log.Base(x, dev)
    }
  }
)

#' Plot histograms of the parameters and weights
#' 
#' @description 
#' Calling [plotHist()] plots histograms of the sampled parameters and weights 
#' from MCMC sampling.
#' 
#' Note, this method is so far only implemented for mictures of Poisson and 
#' Binomial distributions.
#' 
#' @param x An `mcmcoutputpermbase` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Histograms of the MCMC samples.
#' @exportMethod plotHist
#' @describeIn mcmcoutputpermbase-class Plot histograms of the parameters and 
#'   weights
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
#' # Do not use a hierarchical prior.
#' setHier(f_prior) <- FALSE
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' f_outputperm <- mcmcpermute(f_output)
#' plotHist(f_outputperm)
#' }
#' 
#' @seealso 
#' * [mixturemcmc()] for performing MCMC sampling
#' * [mcmcpermute()] for permuting MCMC samples
#' * [plotTraces()] for plotting the traces of sampled values
#' * [plotDens()] for plotting densities of sampled values
#' * [plotSampRep()] for plotting sampling representations of sampled values
#' * [plotPointProc()] for plotting point processes for sampled values
#' * [plotPostDens()] for plotting the posterior density of component parameters
setMethod(
  "plotHist", signature(
    x = "mcmcoutputpermbase",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    dist <- x@model@dist
    if (dist == "poisson") {
      .permhist.Poisson.Base(x, dev)
    } else if (dist == "binomial") {
      .permhist.Binomial.Base(x, dev)
    }
  }
)

#' Plot densities of the parameters and weights
#' 
#' @description 
#' Calling [plotDens()] plots densities of the sampled parameters and weights 
#' from MCMC sampling. 
#' 
#' @param x An `mcmcoutputpermbase` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Densities of the MCMC samples.
#' @exportMethod plotDens
#' @describeIn mcmcoutputpermbase-class
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
#' # Do not use a hierarchical prior.
#' setHier(f_prior) <- FALSE
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' f_outputperm <- mcmcpermute(f_output)
#' plotDens(f_outputperm)
#' 
#' @seealso 
#' * [mixturemcmc()] for performing MCMC sampling
#' * [mcmcpermute()] for permuting MCMC samples
#' * [plotTraces()] for plotting the traces of sampled values
#' * [plotHist()] for plotting histograms of sampled values
#' * [plotSampRep()] for plotting sampling representations of sampled values
#' * [plotPointProc()] for plotting point processes for sampled values
#' * [plotPostDens()] for plotting the posterior density of component parameters
setMethod(
  "plotDens", signature(
    x = "mcmcoutputpermbase",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    dist <- x@model@dist
    if (dist == "poisson") {
      .permdens.Poisson.Base(x, dev)
    } else if (dist == "binomial") {
      .permdens.Binomial.Base(x, dev)
    }
  }
)

#' Plot point processes of the component parameters
#' 
#' @description 
#' Calling [plotPointProc()] plots point processes of the sampled component 
#' parameters from MCMC sampling.  
#' 
#' Note, this method is only implemented for mixtures of Poisson and Binomial 
#' distributions.
#' 
#' @param x An `mcmcoutputpermbase` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Point process of the MCMC samples.
#' @exportMethod plotPointProc
#' @describeIn mcmcoutputpermbase-class Plots point process for the component 
#'   parameters
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
#' # Do not use a hierarchical prior.
#' setHier(f_prior) <- FALSE
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' f_outputperm <- mcmcpermute(f_output)
#' plotPointProc(f_outputperm)
#' }
#' 
#' @seealso 
#' * [mixturemcmc()] for performing MCMC sampling
#' * [mcmcpermute()] for permuting MCMC samples
#' * [plotTraces()] for plotting the traces of sampled values
#' * [plotHist()] for plotting histograms of sampled values
#' * [plotDens()] for plotting densities of sampled values
#' * [plotSampRep()] for plotting sampling representations of sampled values
#' * [plotPostDens()] for plotting posterior densities for sampled values
setMethod(
  "plotPointProc", signature(
    x = "mcmcoutputpermbase",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    dist <- x@model@dist
    if (dist == "poisson") {
      .permpointproc.Poisson(x, dev)
    } else if (dist == "binomial") {
      .permpointproc.Binomial(x, dev)
    }
  }
)

#' Plot sampling representations for the component parameters
#' 
#' @description 
#' Calling [plotSampRep()] plots sampling representations of the sampled 
#' component parameters from MCMC sampling.  
#' 
#' Note, this method is only implemented for mixtures of Poisson and Binomial 
#' distributions.
#' 
#' @param x An `mcmcoutputpermbase` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Sampling representation of the MCMC samples.
#' @exportMethod plotSampRep
#' @describeIn mcmcoutputpermbase-class Plots sampling representations of the 
#'   component parameters
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
#' # Do not use a hierarchical prior.
#' setHier(f_prior) <- FALSE
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' f_outputperm <- mcmcpermute(f_output)
#' plotSampRep(f_outputperm)
#' }
#' 
#' @seealso 
#' * [mixturemcmc()] for performing MCMC sampling
#' * [mcmcpermute()] for permuting MCMC samples
#' * [plotTraces()] for plotting the traces of sampled values
#' * [plotHist()] for plotting histograms of sampled values
#' * [plotDens()] for plotting densities of sampled values
#' * [plotPointProc()] for plotting point processes of sampled values
#' * [plotPostDens()] for plotting posterior densities for sampled values
setMethod(
  "plotSampRep", signature(
    x = "mcmcoutputpermbase",
    dev = "ANY"
  ),
  function(x, dev, ...) {
    dist <- x@model@dist
    if (dist == "poisson") {
      .permsamprep.Poisson(x, dev)
    } else if (dist == "binomial") {
      .permsamprep.Binomial(x, dev)
    }
  }
)

#' Plot posterior densities of the component parameters
#' 
#' @description 
#' Calling [plotPostDens()] plots posterior densities of the sampled component 
#' parameters from MCMC sampling, if the number of components is two. 
#' 
#' Note, this method is so far only implemented for Poisson or Binomial 
#' mixture distributions.
#' 
#' @param x An `mcmcoutputpermbase` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Posterior densities of the MCMC samples.
#' @exportMethod plotPostDens
#' @describeIn mcmcoutputpermbase-class Plots the posterior density of the 
#'   component parameters
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
#' # Do not use a hierarchical prior.
#' setHier(f_prior) <- FALSE
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' f_outputperm <- mcmcpermute(f_output)
#' plotPostDens(f_outputperm)
#' }
#' 
#' @seealso 
#' * [mixturemcmc()] for performing MCMC sampling
#' * [mcmcpermute()] for permuting MCMC samples
#' * [plotTraces()] for plotting the traces of sampled values
#' * [plotHist()] for plotting histograms of sampled values
#' * [plotDens()] for plotting densities of sampled values
#' * [plotSampRep()] for plotting sampling representations of sampled values
#' * [plotPointProc()] for plotting point processes for sampled values
setMethod(
  "plotPostDens", signature(
    x = "mcmcoutputpermbase",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    dist <- x@model@dist
    if (dist == "poisson") {
      .permpostdens.Poisson(x, dev)
    } else if (dist == "binomial") {
      .permpostdens.Binomial(x, dev)
    }
  }
)

### Private functions.
### These functions are not exported.

### Plot
### Traces

#' Plots traces of Poisson mixture samples
#' 
#' @description 
#' For internal usage only. This function plots the traces for sampled values 
#' from a Poisson mixture model.
#' 
#' @param x An `mcmcoutputpermbase` object containing all samples.
#' @param dev A logical indicating if the plot should be shown by a graphical 
#'   device.
#' @return A plot of the traces of sampled values.
#' @noRd
#' 
#' @seealso 
#' * [plotTraces()] for the calling function
".permtraces.Poisson.Base" <- function(x, dev) {
  K <- x@model@K
  trace.n <- K * 2 - 1
  if (.check.grDevice() && dev) {
    dev.new(title = "Traceplots (permuted)")
  }
  par(
    mfrow = c(trace.n, 1), mar = c(1, 0, 0, 0),
    oma = c(4, 5, 4, 4)
  )
  lambda <- x@parperm$lambda
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
  weight <- x@weightperm
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
  axis(1)
  mtext(side = 1, "Iterations", cex = 0.7, line = 3)
}

#' Plots traces of Binomial mixture samples
#' 
#' @description 
#' For internal usage only. This function plots the traces for sampled values 
#' from a Binomial mixture model.
#' 
#' @param x An `mcmcoutputpermbase` object containing all samples.
#' @param dev A logical indicating if the plot should be shown by a grapical 
#'   device.
#' @return A plot of the traces of sampled values.
#' @noRd
#' 
#' @seealso 
#' * [plotTraces()] for the calling function
".permtraces.Binomial.Base" <- function(x, dev) {
  K <- x@model@K
  trace.n <- K * 2 - 1
  if (.check.grDevice() && dev) {
    dev.new(title = "Traceplots (permuted)")
  }
  par(
    mfrow = c(trace.n, 1), mar = c(1, 0, 0, 0),
    oma = c(4, 5, 4, 4)
  )
  p <- x@parperm$p
  for (k in 1:K) {
    plot(p[, k],
      type = "l", axes = F,
      col = "gray20", xlab = "", ylab = ""
    )
    axis(2, las = 2, cex.axis = 0.7)
    mtext(
      side = 2, las = 2, bquote(p[k = .(k)]),
      cex = 0.6, line = 3
    )
  }
  weight <- x@weightperm
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
  axis(1)
  mtext(side = 1, "Iterations", cex = 0.7, line = 3)
}

#' Plots traces of exponential mixture samples
#' 
#' @description 
#' For internal usage only. This function plots the traces for sampled values 
#' from a exponential mixture model.
#' 
#' @param x An `mcmcoutputpermbase` object containing all samples.
#' @param dev A logical indicating if the plot should be shown by a graphical 
#'   device.
#' @return A plot of the traces of sampled values.
#' @noRd
#' 
#' @seealso 
#' * [plotTraces()] for the calling function
".permtraces.Exponential.Base" <- function(x, dev) {
  K <- x@model@K
  trace.n <- K * 2 - 1
  if (.check.grDevice() && dev) {
    dev.new(title = "Traceplots")
  }
  par(
    mfrow = c(trace.n, 1), mar = c(1, 0, 0, 0),
    oma = c(4, 5, 4, 4)
  )
  lambda <- x@parperm$lambda
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

#' Plots traces of sampled weights
#' 
#' @description 
#' For internal usage only. This function plots the traces for sampled values 
#' of the weights for any mixture model.
#' 
#' @param x An `mcmcoutputpermbase` object containing all samples.
#' @param dev A logical indicating if the plot should be shown by a grapical 
#'   device.
#' @return A plot of the traces of sampled values.
#' @noRd
#' 
#' @seealso 
#' * [plotTraces()] for the calling function
".permtraces.Weights.Base" <- function(x, dev, col) {
  weight <- x@weightperm
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

### Traces log-likelihoods: Plots traces for the log-likelihoods.
#' Plots traces of log-likelihood samples
#' 
#' @description 
#' For internal usage only. This function plots the traces for the 
#' log-likelihoods of sampled values from any mixture model.
#' 
#' @param x An `mcmcoutputpermfix` object containing all samples.
#' @param dev A logical indicating if the plot should be shown by a graphical 
#'   device.
#' @return A plot of the traces of sampled values.
#' @noRd 
#' 
#' @seealso 
#' * [plotTraces()] for the calling function
".permtraces.Log.Base" <- function(x, dev) {
  if (.check.grDevice() && dev) {
    dev.new(title = "Log Likelihood Traceplots (permuted)")
  }
  if (col) {
    cscale <- rainbow(3, start = 0, end = .5)
  } else {
    cscale <- gray.colors(3, start = 0, end = .15)
  }
  par(
    mfrow = c(3, 1), mar = c(1, 0, 0, 0),
    oma = c(4, 5, 4, 4)
  )
  mixlik <- x@logperm$mixlik
  plot(mixlik,
    type = "l", axes = F,
    col = cscale[3], xlab = "", ylab = ""
  )
  axis(2, las = 2, cex.axis = 0.7)
  mtext(
    side = 2, las = 3, "mixlik", cex = 0.6,
    line = 3
  )
  mixprior <- x@logperm$mixprior
  plot(mixprior,
    type = "l", axes = F,
    col = cscale[2], xlab = "", ylab = ""
  )
  axis(2, las = 2, cex.axis = 0.7)
  mtext(
    side = 2, las = 3, "mixprior", cex = 0.6,
    line = 3
  )
  cdpost <- x@logperm$cdpost
  plot(mixprior,
    type = "l", axes = F,
    col = cscale[3], xlab = "", ylab = ""
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
#' @param x An `mcmcoutputpermbase` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with histograms for the sampled parameters and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotHist()] for the calling function
".permhist.Poisson.Base" <- function(x, dev) {
  K <- x@model@K
  if (.check.grDevice() && dev) {
    dev.new(title = "Histograms (permuted)")
  }
  lambda <- x@parperm$lambda
  weight <- x@weightperm
  vars <- cbind(lambda, weight[, seq(1:(K - 1))])
  lab.names <- vector("list", 2 * K - 1)
  for (k in 1:K) {
    lab.names[[k]] <- bquote(lambda[.(k)])
  }
  for (k in (K + 1):(2 * K - 1)) {
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
#' @param x An `mcmcoutputpermbase` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with histograms for the sampled parameters and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotHist()] for the calling function
".permhist.Binomial.Base" <- function(x, dev) {
  K <- x@model@K
  if (.check.grDevice() && dev) {
    dev.new(title = "Histograms (permuted)")
  }
  p <- x@parperm$p
  weight <- x@weightperm
  vars <- cbind(p, weight[, seq(1:(K - 1))])
  lab.names <- vector("list", 2 * K - 1)
  for (k in 1:K) {
    lab.names[[k]] <- bquote(p[.(k)])
  }
  for (k in (K + 1):(2 * K - 1)) {
    lab.names[[k]] <- bquote(eta[.(k - K)])
  }
  .symmetric.Hist(vars, lab.names)
}

### Densities

#' Plot densities of Poisson samples
#' 
#' @description 
#' For internal usage only. This function plots densities of sampled Poisson 
#' parameters and weights.
#' 
#' @param x An `mcmcoutputpermbase` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with densities for the sampled parameters and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotDens()] for the calling function
".permdens.Poisson.Base" <- function(x, dev) {
  K <- x@model@K
  if (.check.grDevice() && dev) {
    dev.new(title = "Histograms (permuted)")
  }
  lambda <- x@parperm$lambda
  weight <- x@weightperm
  vars <- cbind(lambda, weight[, seq(1:(K - 1))])
  lab.names <- vector("list", 2 * K - 1)
  for (k in 1:K) {
    lab.names[[k]] <- bquote(lambda[.(k)])
  }
  for (k in (K + 1):(2 * K - 1)) {
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
#' @param x An `mcmcoutputpermbase` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with densities for the sampled parameters and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotDens()] for the calling function
".permdens.Binomial.Base" <- function(x, dev) {
  K <- x@model@K
  if (.check.grDevice() && dev) {
    dev.new(title = "Histograms (permuted)")
  }
  p <- x@parperm$p
  weight <- x@weightperm
  vars <- cbind(p, weight[, seq(1:(K - 1))])
  lab.names <- vector("list", 2 * K - 1)
  for (k in 1:K) {
    lab.names[[k]] <- bquote(p[.(k)])
  }
  for (k in (K + 1):(2 * K - 1)) {
    lab.names[[k]] <- bquote(eta[.(k - K)])
  }
  .symmetric.Dens(vars, lab.names)
}