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
#' @param An `array` of dimension `N x storeS` containing the last 
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
#' The mcmcoutputperm class stores MCMC samples after relabeling (permuting). 
#' 
#' @details 
#' Calling [mcmcpermute()] on an `mcmcoutput` class permutes the labels of the 
#' components and generates an object of class `mcmcoutputperm`. Note, the 
#' number of samples of the `mcmcoutputperm` object could be less than the 
#' original number of MCMC samples due to some samples where both components 
#' get assigned to the same label and henceforth get eliminated from further 
#' analysis.
#' 
#' This class union includes all classes that define objects for permuted 
#' MCMC samples and is used to dispatch methods for `mcmcoutputperm` objects.
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