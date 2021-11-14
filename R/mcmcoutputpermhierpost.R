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

#' Finmix `mcmcoutputpermhierpost` class 
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
#' @exportClass mcmcoutputpermhierpost
#' @rdname mcmcoutputpermhierpost-class
#' @keywords internal
#' 
#' @seealso 
#' * [mcmcoutputbase-class] for the parent class
#' * [mcmcpermind-class] for the parent class
#' * [mcmcpermute()] for performing permutation of MCMC samples
.mcmcoutputpermhierpost <- setClass("mcmcoutputpermhierpost",
  contains = c(
    "mcmcpermindhier",
    "mcmcpermindpost",
    "mcmcoutputhierpost"
  ),
  validity = function(object) {
    ## else: OK
    TRUE
  }
)

#' Initializer of the `mcmcoutputpermhierpost` class
#' 
#' @description
#' Only used implicitly. The initializer stores the data into the slots of the 
#' passed-in object.
#' 
#' @param .Object An object: see the "initialize Methods" section in 
#'   [initialize].
#' @param mcmcoutput An `mcmcoutput` class containing the results from MCMC 
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
#' @param hyperperm A named list containing the (permuted) parameters of the 
#'   hierarchical prior.
#' @param postperm A named list containing a named list `par` with array(s) of 
#'   parameters from the posterior density. 
#' @param entropyperm An `array` of dimension `Mperm x 1` containing the 
#'   entropy for each MCMC permuted draw.
#' @param STperm An `array` of dimension `Mperm x 1` containing all permuted 
#'   MCMC states, for the last observation in slot `@@y` of the `fdata` object 
#'   passed in to [mixturemcmc()] where a state is defined for non-Markov 
#'   models as the last indicator of this observation.  
#' @param Sperm An `array` of dimension `N x storeS` containing the last 
#'   `storeS` permuted indicators. `storeS` is defined in the slot `@@storeS` 
#'   of the `mcmc` object passed into [mixturemcmc()].
#' @param NKperm An `array` of dimension `Mperm x K` containing the numbers 
#'   of observations assigned to each component.
#' 
#' @keywords internal
#' 
#' @seealso 
#' * [Classes_Details] for details of class definitions, and 
#' * [setOldClass] for the relation to S3 classes
setMethod(
  "initialize", "mcmcoutputpermhierpost",
  function(.Object, mcmcoutput, Mperm = integer(),
           parperm = list(), relabel = character(),
           weightperm = array(), logperm = list(),
           hyperperm = list(), postperm = list(), 
           entropyperm = array(), STperm = array(), 
           Sperm = array(), NKperm = array()) {
    .Object@M <- mcmcoutput@M
    .Object@burnin <- mcmcoutput@burnin
    .Object@ranperm <- mcmcoutput@ranperm
    .Object@par <- mcmcoutput@par
    .Object@weight <- mcmcoutput@weight
    .Object@log <- mcmcoutput@log
    .Object@hyper <- mcmcoutput@hyper
    .Object@post <- mcmcoutput@post
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
    .Object@hyperperm <- hyperperm
    .Object@postperm <- postperm
    .Object@entropyperm <- entropyperm
    .Object@STperm <- STperm
    .Object@Sperm <- Sperm
    .Object@NKperm <- NKperm
    .Object
  }
)

#' Shows a summary of an `mcmcoutputpermhierpost` object.
#' 
#' @description
#' Calling [show()] on an `mcmcoutputpermhierpost` object gives an overview 
#' of the `mcmcoutputpermhierpost` object.
#' 
#' @param object An `mcmcoutputpermhierpost` object.
#' @returns A console output listing the slots and summary information about
#'   each of them. 
#' @exportMethod show
#' @keywords internal
setMethod(
  "show", "mcmcoutputpermhierpost",
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
      "     weight      :",
      paste(dim(object@weight), collapse = "x"), "\n"
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
      "     post        : List of",
      length(object@post), "\n"
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
      "     hyperperm   : List of",
      length(object@hyperperm), "\n"
    )
    cat(
      "     postperm    : List of",
      length(object@postperm), "\n"
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
#' If `lik` is set to `0` the parameters of the components, the posterior 
#' parameters, and the parameters of the hierarchical prior are plotted 
#' together with `K-1` weights.
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
#' @exportMethod  plotTraces
#' @keywords internal
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
#' plotTraces(f_outputperm, lik = 0)
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
    x = "mcmcoutputpermhierpost",
    dev = "ANY",
    lik = "ANY",
    col = "ANY"
  ),
  function(x, dev = TRUE, lik = 1, col = FALSE, ...) {
    dist <- x@model@dist
    if (lik %in% c(0, 1)) {
      if (dist == "poisson") {
        .permtraces.Poisson.Base.Hier(x, dev)
      } else if (dist == "binomial") {
        .permtraces.Binomial.Base(x, dev)
      } else if (dist == "exponential") {
        .permtraces.Exponential.Base(x, dev)
        .permtraces.Weights.Base(x, dev, col)
      } else if (dist == "normal") {
        .permtraces.Normal.Hier(x, dev)
        .permtraces.Weights.Base(x, dev, col)
      } else if (dist == "student") {
        .permtraces.Student.Hier(x, dev)
        .permtraces.Weights.Base(x, dev, col)
      } else if (dist == "normult") {
        .permtraces.Normult.Hier(x, dev, col)
        .permtraces.Weights.Base(x, dev, col)
      } else if (dist == "studmult") {
        .permtraces.Studmult.Hier(x, dev, col)
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
#' from MCMC sampling. In addition the parameters of the hierarchical prior are 
#' plotted.
#' 
#' Note, this method is so far only implemented for mictures of Poisson and 
#' Binomial distributions.
#' 
#' @param x An `mcmcoutputpermhierpost` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Histograms of the MCMC samples.
#' @exportMethod plotHist
#' @keywords internal
#' 
#' @examples 
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc()
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' f_outputperm <- mcmcpermute(f_output)
#' plotHist(f_outputperm)
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
    x = "mcmcoutputpermhierpost",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    dist <- x@model@dist
    if (dist == "poisson") {
      .permhist.Poisson.Base.Hier(x, dev)
    } else if (dist == "binomial") {
      .permhist.Binomial.Base(x, dev)
    }
  }
)

#' Plot densities of the parameters and weights
#' 
#' @description 
#' Calling [plotDens()] plots densities of the sampled parameters and weights 
#' from MCMC sampling. In addition, the parameters of the hierarchical prior 
#' are plotted.
#' 
#' Note, this method is so far only implemented for mixtures of Poisson and 
#' Binomial distributions. 
#' 
#' @param x An `mcmcoutputpermhierpost` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Densities of the MCMC samples.
#' @exportMethod plotDens
#' @keywords internal
#' 
#' @examples 
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc()
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
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
    x = "mcmcoutputpermhierpost",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    dist <- x@model@dist
    if (dist == "poisson") {
      .permdens.Poisson.Base.Hier(x, dev)
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
#' @param x An `mcmcoutputpermhierpost` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Densities of the MCMC samples.
#' @exportMethod plotPointProc
#' @keywords internal
#' 
#' @examples 
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc()
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' f_outputperm <- mcmcpermute(f_output)
#' plotPointProc(f_outputperm)
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
    x = "mcmcoutputpermhierpost",
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
#' @param x An `mcmcoutputpermhierpost` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Densities of the MCMC samples.
#' @exportMethod plotSampRep
#' @keywords internal
#' 
#' @examples 
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc()
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' f_outputperm <- mcmcpermute(f_output)
#' plotSampRep(f_outputperm)
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
    x = "mcmcoutputpermhierpost",
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
#' @param x An `mcmcoutputpermhierpost` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Posterior densities of the MCMC samples.
#' @exportMethod plotPostDens
#' @keywords internal
#' 
#' @examples 
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc()
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' f_outputperm <- mcmcpermute(f_output)
#' plotPostDens(f_outputperm)
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
    x = "mcmcoutputpermhierpost",
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

#' Finmix `mcmcoutputperm` class
#' 
#' @description 
#' The `mcmcoutputperm` class stores MCMC samples after relabeling (permuting). 
#' 
#' @details 
#' Calling [mcmcpermute()] on an `mcmcoutput` class permutes the labels of the 
#' components and generates an object of class `mcmcoutputperm`. Note, the 
#' number of samples of the `mcmcoutputperm` object could be less than the 
#' original number of MCMC samples due to some samples where both components 
#' get assigned to the same label and henceforth get eliminated from further 
#' analysis.
#' 
#' The class `mcmcoutputperm` is a class union and includes all classes that 
#' define objects for permuted MCMC samples and is used to dispatch methods for 
#' `mcmcoutputperm` objects. For the user this detail is not important, 
#' especially as this class has no exported constructor. Objects are solely 
#' constructed internally within the function [mcmcpermute()]. 
#' 
#' An object of class `mcmcoutputperm` inherits all slots from its parent class 
#' [mcmcoutput][mcmcoutput-class]. In addition it contains slots that store 
#' data from permutation. These slots are listed below
#' 
#' ## Class methods
#' Similar to the contained classes [mcmcoutput][mcmcoutput-class] this class
#' comes along with a couple of methods that should give the user some comfort
#' in handling the permuted sampling results. There are no setters for this
#' class as the slots are only set internally.
#' 
#' ### Show
#' * `show()` shows a short summary of the object's slots.
#' 
#' ### Getters
#' * `getMperm()` returns the `Mperm` slot.
#' * `getParperm()` returns the `parperm` slot. 
#' * `getLogperm()` returns the `parperm` slot. 
#' * `getHyperperm()` returns the `hyperparm` slot. 
#' * `getPostperm()` returns the `postperm` slot. 
#' * `getEntropyperm()` returns the `entropyperm` slot.
#' * `getSTperm()` returns the `STperm` slot.
#' * `getSperm()` returns the `Sperm` slot.
#' * `getNKperm()` returns the `NKperm` slot. 
#' 
#' ### Plotting
#' Plotting functionality for the `mcmcoutputperm` class is so far only 
#' implemented for mixtures of Binomial or Poisson distributions.
#' 
#' * `plotTraces()` plots traces of relabeled MCMC sampling. See [plotTraces()] 
#'   for further information.
#' * `plotHist()` plots histograms of relabeled parameters and weights. See 
#'   [plotHist()] for further information.
#' * `plotDens()` plots densities of relabeled parameters and weights. See 
#'   [plotDens()] for further information.
#' * `plotPointProc()` plots the point process of relabeled component 
#'   parameters. See [plotPointProc] for further information.
#' * `plotSampRep()` plots the sampling representation of relabeled component 
#'   parameters. See [plotSampRep()] for further information. 
#' * `plotPostDens()` plots the posterior density of component parameters. Note 
#'   that this function can only be applied for mixtures of two components. See 
#'   [plotPostDens()] for further information.    
#' 
#' * `Mperm` An integer defining the number of permuted MCMC samples.
#' * `parperm` A named list containing the permuted component parameter 
#'   samples from MCMC sampling.
#' * `relabel` A character specifying the relabeling algorithm used for 
#'   permuting the MCMC samples.
#' * `weightperm` An array of dimension `MpermxK` containing the 
#'   relabeled weight parameters. This slot is not available for models with 
#'   fixed indicators as weights do not get sampled for such models.
#' * `logperm` A named list containing the mixture log-likelihood, the 
#'   prior log-likelihood, and for models with unknown indicators the complete 
#'   data posterior log-likelihood for the permuted MCMC samples.
#' * `hyperperm` A named list containing the (permuted) parameters of the 
#'   hierarchical prior. This slot is only available, if a hierarchical prior 
#'   had been used for sampling, i.e. the slot `hier` of the 
#'   [prior][prior-class] had been set to `TRUE`.
#' * `postperm` A named list containing a named list `par` with array(s) of 
#'   parameters from the posterior density. This slot is only available if 
#'   the hyperparameter `storepost` in the [mcmc][mcmc-class] object had been 
#'   set to `TRUE`.
#' * `entropyperm` An `array` of dimension `Mpermx1` containing the 
#'   entropy for each MCMC permuted draw. This slot is only available for 
#'   models with unknown indicators.
#' `STperm` An `array` of dimension `Mpermx1` containing all permuted 
#'   MCMC states, for the last observation in slot `y` of the `fdata` object 
#'   passed in to [mixturemcmc()] where a state is defined for non-Markov 
#'   models as the last indicator of this observation. This slot is only 
#'   available for models with unknown indicators.  
#' * `Sperm` An `array` of dimension `N x storeS` containing the last 
#'   `storeS` permuted indicators. `storeS` is defined in the slot `storeS` 
#'   of the `mcmc` object passed into [mixturemcmc()]. This slot is only 
#'   available for models with unknown indicators.
#' * `NKperm` An `array` of dimension `Mperm x K` containing the numbers 
#'   of observations assigned to each component. This slot is only available for 
#'   models with unknown indicators.
#'  
#' @exportClass mcmcoutputperm
#' @rdname mcmcoutputperm-class
setClassUnion(
  "mcmcoutputperm",
  c(
    "mcmcoutputpermfix",
    "mcmcoutputpermfixhier",
    "mcmcoutputpermfixpost",
    "mcmcoutputpermfixhierpost",
    "mcmcoutputpermbase",
    "mcmcoutputpermhier",
    "mcmcoutputpermpost",
    "mcmcoutputpermhierpost"
  )
)

#' Plots traces of MCMC sampling
#' 
#' @description 
#' `plotTraces()` is a class method for [mcmcoutput][mcmcoutput-class] and 
#' [mcmcoutputperm][mcmcoutputperm-class] objects. For the former class it 
#' plots the traces of MCMC samples and for the latter of the corresponding 
#' permuted samples coming from relabeling. 
#' 
#' @details
#' Calling [plotTraces()] with `lik` set to `1`, plots the MCMC traces of the 
#' mixture log-likelihood, the mixture log-likelihood of the prior 
#' distribution, or the log-likelihood of the complete data posterior, if the 
#' model has unknown indicators. 
#' 
#' If `lik` is set to `0` the parameters of the components, the posterior 
#' parameters, and the parameters of the hierarchical prior are plotted 
#' together with `K-1` weights.
#' 
#' ## Hierarchical priors
#' In case of hierarchical priors, the function also plots traces from the 
#' sampled hierarchical prior's parameters, in case `lik` is set to `1`.
#' 
#' ## Posterior density parameters
#' In case posterior density parameters had been stored in MCMC sampling, the 
#' traces of these parameters are added to the plot. 
#' 
#' @param x An `mcmcoutput` or `mcmcoutputperm` object containing all sampled 
#'   values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param lik An integer indicating, if the log-likelihood traces should be 
#'   plotted (default). If set to `0` the traces for the parameters 
#'   and weights are plotted instead. 
#' @param col A logical indicating, if the plot should be colored.
#' @param ... Further arguments to be passed to the plotting function.
#' @return A plot of the traces of the MCMC samples.
#' @name plotTraces
#' @rdname plotTraces-method
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
#' plotTraces(f_outputperm, lik = 0)
#' 
#' @seealso 
#' * [mixturemcmc()] for performing MCMC sampling
#' * [mcmcpermute()] for permuting MCMC samples 
#' * [plotHist()] for plotting histograms of sampled values
#' * [plotDens()] for plotting densities of sampled values
#' * [plotSampRep()] for plotting sampling representations of sampled values
#' * [plotPointProc()] for plotting point processes for sampled values
#' * [plotPostDens()] for plotting the posterior density of component parameters
#' * [mcmcoutput-class] for the class definition of `mcmcoutput`
#' * [mcmcoutputperm-class] for the class definition of `mcmcoutputperm` 
NULL

#' Plot histograms of the parameters and weights
#' 
#' @description 
#' `plotHist()` is a class method for [mcmcoutput][mcmcoutput-class] and 
#' [mcmcoutputperm][mcmcoutputperm-class] objects. For the former class it 
#' plots histograms of MCMC samples and for the latter of the corresponding 
#' permuted samples coming from relabeling. 
#' 
#' @details
#' Calling [plotHist()] plots histograms of the sampled parameters and weights 
#' from MCMC sampling. Note, for relabeled MCMC samples this method is so far 
#' only implemented for mixtures of Poisson and Binomial distributions.
#' 
#' ## Hierarchical priors
#' In case that hierarchical priors had been used in MCMC sampling histograms 
#' of the sampled parameters of the hierarchical prior are added to the plot.
#' 
#' ## Posterior density parameters
#' In case that posterior density parameters had been stored in MCMC sampling, 
#' histograms of these parameters are added to the plot. 
#' 
#' @param x An `mcmcoutput` or `mcmcoutputperm` object containing all sampled 
#'   values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Histograms of the MCMC samples.
#' @name plotHist
#' @rdname plotHist-method
#' 
#' @examples 
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc()
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' f_outputperm <- mcmcpermute(f_output)
#' plotHist(f_outputperm)
#' 
#' @seealso 
#' * [mixturemcmc()] for performing MCMC sampling
#' * [mcmcpermute()] for permuting MCMC samples
#' * [plotTraces()] for plotting the traces of sampled values
#' * [plotDens()] for plotting densities of sampled values
#' * [plotSampRep()] for plotting sampling representations of sampled values
#' * [plotPointProc()] for plotting point processes for sampled values
#' * [plotPostDens()] for plotting the posterior density of component parameters
#' * [mcmcoutput-class] for the class definition of `mcmcoutput`
#' * [mcmcoutputperm-class] for the class definition of `mcmcoutputperm` 
NULL

#' @title Plot densities of the parameters and weights
#' 
#' @description 
#' `plotDens()` is a class method for [mcmcoutput][mcmcoutput-class] and 
#' [mcmcoutputperm][mcmcoutputperm-class] objects. For the former class it 
#' plots densities of MCMC samples and for the latter of the corresponding 
#' permuted samples coming from relabeling. 
#' 
#' @details
#' Calling [plotDens()] plots densities of the sampled parameters and weights 
#' from MCMC sampling. Note, for relabeled MCMC samples this method is so far 
#' only implemented for mixtures of Poisson and Binomial distributions.
#' 
#' ## Hierarchical priors
#' In case that hierarchical priors had been used in MCMC sampling densities 
#' of the sampled parameters of the hierarchical prior are added to the plot.
#' 
#' ## Posterior density parameters
#' In case that posterior density parameters had been stored in MCMC sampling, 
#' densities of these parameters are added to the plot. 
#' 
#' @param x An `mcmcoutput` or `mcmcoutputperm` object containing all sampled 
#'   values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Densities of the MCMC samples.
#' @name plotDens
#' @rdname plotDens-method
#' 
#' @examples 
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc()
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' f_outputperm <- mcmcpermute(f_output)
#' plotDens(f_outputperm)
#' 
#' @seealso 
#' * [mixturemcmc()] for performing MCMC sampling
#' * [mcmcpermute()] for permuting MCMC samples
#' * [plotTraces()] for plotting the traces of sampled values
#' * [plotDens()] for plotting densities of sampled values
#' * [plotSampRep()] for plotting sampling representations of sampled values
#' * [plotPointProc()] for plotting point processes for sampled values
#' * [plotPostDens()] for plotting the posterior density of component parameters
#' * [mcmcoutput-class] for the class definition of `mcmcoutput`
#' * [mcmcoutputperm-class] for the class definition of `mcmcoutputperm` 
NULL

#' Plot the point process of the component parameters
#' 
#' @description 
#' Calling [plotPointProc()] on an object of class `mcmcoutput` or 
#' `mcmcoutputperm` plots the point process of the sampled component parameters 
#' from MCMC sampling, either the original parameters or the relabeled ones.
#' 
#' @details 
#' The point process is used to identify the number of components in the 
#' underlying distribution of the data for mixtures with unknown number of 
#' components (see Frühwirth-Schnatter (2006, Subsection 3.7.1)). The number of 
#' clusters that evolve in the plot give a hint on the true number of 
#' components in the mixture distribution. The MCMC draws will scatter around 
#' the points corresponding to the true point process of the mixture model. The 
#' spread of the clusters represent the uncertainty of estimating the points. 
#' 
#' For mixtures with univariate component parameters (e.g. Poisson, 
#' Exponential) the component parameters are plotted against draws from a 
#' standard normal distribution. For mixtures with bivariate component 
#' parameters (e.g. Normal) the first parameters are plotted against the 
#' second ones. For mixtures with multivariate component parameters a point 
#' process for each type of mixture model is plotted. 
#' 
#' Note that this method for `mcmcoutputperm` objects is only implemented for 
#' mixtures of Poisson and Binomial distributions.
#' 
#' @param x An `mcmcoutput` or `mcmcoutputperm` object containing all sampled 
#'   values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return The point process of the MCMC samples.
#' @rdname plotPointProc-method
#' @name plotPointProc
#' 
#' @examples 
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc()
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' f_outputperm <- mcmcpermute(f_output)
#' plotPointProc(f_outputperm)
#' 
#' @seealso 
#' * [mixturemcmc()] for performing MCMC sampling
#' * [mcmcpermute()] for permuting MCMC samples
#' * [plotTraces()] for plotting the traces of sampled values
#' * [plotHist()] for plotting histograms of sampled values
#' * [plotDens()] for plotting densities of sampled values
#' * [plotSampRep()] for plotting sampling representations of sampled values
#' * [plotPostDens()] for plotting posterior densities for sampled values
NULL

#' Plot the sampling representation of component parameters
#' 
#' @description
#' Calling [plotSampRep()] on an object of class `mcmcoutput` or 
#' `mcmcoutputperm` plots the sampling representation of the sampled component 
#' parameters from MCMC sampling, either the original parameters or the 
#' relabeled ones (`mcmcoutputperm`).
#' 
#' @details 
#' To visualize the posterior density of the component parameters the MCMC 
#' draws are used as a sampling representation. Each combination of component 
#' parameters is plotted in a scatter to visualize the contours of the 
#' posterior density. For bivariate component parameters this could also be 
#' done by estimating and plotting the density directly, but for 
#' higher-dimensional parameter vectors this is not anymore possible and so 
#' sampling representations define a proper solution for visualization and 
#' allow us to study how a specific dimension of the parameter vector differs 
#' among the various components of the mixture distribution. If this element 
#' is significantly different among components we will observe `K(K-1)` modes 
#' in the sampling representation. On the other side, if this element is 
#' mainly the same among the components of the mixture, we will rather observe 
#' a single cluster.
#' 
#' As Frühwirth-Schnatter (2006) writes, "One informal method for diagnosing 
#' mixtures is mode hunting in the mixture posterior density 
#' (Frühwirth-Schnatter, 2001b). It is based on the observation that with an 
#' increasing number of observations, the mixture likelihood function has `K!` 
#' dominant modes if the data actually arise from a finite mixture distribution 
#' with `K` components, and that less than `K!` dominant modes are present if 
#' the finite mixture model is overfitting the number of components." The 
#' sampling representation helps to perform this mode hunting in practice.
#' 
#' Note that this method for `mcmcoutputperm` objects is only implemented for 
#' mixtures of Poisson and Binomial distributions.
#' 
#' @param x An `mcmcoutput` or `mcmcoutputperm` object containing all sampled 
#'   values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return The sampling representation of the MCMC samples.
#' @rdname plotSampRep-method
#' @name plotSampRep
#' 
#' @examples 
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc()
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' plotSampRep(f_output)
#' 
#' @references 
#' * Frühwirth-Schnatter (2006), "Finite Mixture and Markov Switching Models"
#' * Frühwirth-Schnatter, S. (2001b), "Markov chain Monte Carlo estimation of 
#'   classical and dynamic switching and mixture models." Journal of the 
#'   American Statistical Association 96, 194–209.
#'   
#' @seealso 
#' * [mixturemcmc()] for performing MCMC sampling
#' * [mcmcpermute()] for permuting MCMC samples
#' * [plotTraces()] for plotting the traces of sampled values
#' * [plotHist()] for plotting histograms of sampled values
#' * [plotDens()] for plotting densities of sampled values
#' * [plotPointProc()] for plotting the point process of sampled values
#' * [plotPostDens()] for plotting posterior densities for sampled values
NULL

#' Plot the posterior density of component parameters
#' 
#' @description
#' Calling [plotPostDens()] on an object of class `mcmcoutput` or 
#' `mcmcoutputperm` plots the posterior density of the sampled component 
#' parameters from MCMC sampling, either the original parameters or the 
#' relabeled ones (`mcmcoutputperm`).
#' 
#' @details 
#' Next to sampling representations and the point process of MCMC samples the 
#' posterior density of component parameters can also be plotted directly for 
#' finite mixture distributions with ` K=2` components and a single parameter.
#' The posterior density will always be bimodal due to to label-switching in 
#' the MCMC sampling. This could change when considering a relabeld MCMC sample 
#' (`mcmcoutputperm` object). 
#' 
#' Note that this method for `mcmcoutputperm` objects is only implemented for 
#' mixtures of Poisson and Binomial distributions.
#' 
#' @param x An `mcmcoutput` or `mcmcoutputperm` object containing all sampled 
#'   values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return The posterior density of the MCMC samples.
#' @rdname plotPostDens-method
#' @name plotPostDens
#' 
#' @examples 
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc()
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' plotPostDens(f_output)
#' 
#' @seealso 
#' * [mixturemcmc()] for performing MCMC sampling
#' * [mcmcpermute()] for permuting MCMC samples
#' * [plotTraces()] for plotting the traces of sampled values
#' * [plotHist()] for plotting histograms of sampled values
#' * [plotDens()] for plotting densities of sampled values
#' * [plotPointProc()] for plotting the point process of sampled values
#' * [plotSampRep()] for plotting the sampling representation for sampled values
NULL

#' Extract sub-chains from MCMC samples
#' 
#' @description 
#' Calling [subseq()] on an `mcmcoutput` or `mcmcoutputperm` object creates a 
#' sub-chain defined by the argument `index`. Sub-chains can be used to further 
#' investigate convergence of MCMC sampling. 
#' 
#' @details 
#' Running MCMC sampling should by time result in a roughly stationary sequence 
#' of random draws. If trace plots do not show this stationary pattern MCMC 
#' sampling should be run with a longer burn-in period until the sampling 
#' distribution has converged. Another possibility is to remove the first draws. 
#' Removing the first draws can be achieved by calling `subseq()` on the object 
#' holding the MCMC samples. 
#' In case of autocorrelations in the traces it is also possible to extract 
#' every `t`-th value by setting the `index` argument accordingly. 
#' 
#' @param object An `mcmcoutput` or `mcmcoutputperm` object containing samples 
#'   from MCMC samples.
#' @param index A logical `array` of dimension `Mx1` defining the schema for 
#'   the sub-chain.
#' @return An `mcmcoutput` or `mcmcoutputperm` object containing the 
#'   sub-chained MCMC samples.
#' @rdname subseq-method
#' @name subseq
#' 
#' @examples 
#' # Define a mixture of Poisson distributions.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2)
#' # Simulate data from the model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc(storepost = FALSE)
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' # Define a sub-chain randomly.
#' index <- array(sample(c(FALSE, TRUE), size = getM(f_output), replace = TRUE))
#' # Extract the sub-chain.
#' subseq(f_output, index)
#' 
#' @seealso 
#' * [mcmcoutput-class] for the class storing MCMC samples
#' * [mcmcoutputperm-class] for the corresponding class for re-labeled MCMC 
#'   samples
#' * [plotTraces()] for plotting traces to be used for a convergence analysis
#' * [swapElements()] for swapping elements in MCMC samples
NULL

#' Swap elements of MCMC samples
#' 
#' @description 
#' Calling `swapElements()` on an `mcmcoutput` object 
#' swaps all labels by the schema given in the `index` argument. 
#' 
#' @details 
#' This function is merely a utility function that simplifies relabeling for 
#' users and developers. For relabeling the labels have to be permuted and 
#' depending on the MCMC sampling chosen there could be a lot of different 
#' slots that need to be permuted. `swapElements()` swaps the elements in any 
#' slot that needs to be relabeled. 
#'  
#' @param object An `mcmcoutput` object containing the 
#'   sampled values.
#' @param index An array specifying the extraction of the values.
#' @return An `mcmcoutput` object with swapped elements.
#' @rdname swapElements-method
#' @name swapElements
#' 
#' @examples 
#' \dontrun{
#' # Generate a model of Poisson distributions.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2)
#' # Simulate data from the model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc(storepost = FALSE)
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' index <- matrix(c(1, 2), nrow = getM(f_output) + 1, 
#'                 ncol = 2)[1:getM(f_output),]
#' swapElements(f_output, index)
#' }
#' 
#' @seealso 
#' * [mcmcoutput-class] for the class definition
#' * [subseq()] for generating sub-chains from MCMC samples
#' * [mcmcpermute()] for a calling function
NULL

#' Extracts single samples from a multivariate Normal mixture
#' 
#' @description 
#' Calling [extract()] on an `mcmcoutput` object with a multivariate Normal 
#' mixture model extracts single samples. 
#' 
#' @details 
#' This function simplifies the analysis of multivariate Normal mixtures that
#' come along with matrices instead of vectors for component parameters as it
#' extracts the mean matrix, the variance matrices and in addition the inverted
#' variance matrices with a single call. In additon, it enriches the output
#' object with metadata like the dimension of the data `r`, the number of
#' components `K`, and the distribution (in this case `"normult`).
#' 
#' @param object An `mcmcoutput` or `mcmcoutputperm` object containing the MCMC 
#'   samples.
#' @param index An `integer` specifying the dimension to extract.
#' @return An `mcmcextract` object containing the parameters, weights, and 
#'   metadata of the extracted dimension. 
#' @rdname extract-method
#' @name extract  
#' 
#' @examples 
#' # Generate a multivariate Normal mixture model.
#' means <- matrix(c(1, 2, 2, 4), nrow = 2) 
#' var1 <- matrix(c(1, 0.3, 0.3, 2), nrow=2)
#' var2 <- matrix(c(3, 0.3, 0.3, 6), nrow=2)
#' vars <- array(c(var1,var2), dim = c(2,2,2))
#' f_model <- model(dist='normult', K = 2, r = 2, par = list(mu=means, sigma=vars))
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc(storepost = FALSE)
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' # Extract a single MCMC sample.
#' f_output1 <- extract(f_output, index = 1000)
#' 
#' @seealso 
#' * [mcmcoutput-class] for the definition of the `mcmcoutput` class
#' * [mcmcoutputperm-class] for the definition of the `mcmcoutputperm` class
#' * [mcmcextract-class] for the output class
NULL