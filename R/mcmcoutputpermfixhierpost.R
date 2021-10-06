# Copyright (C) 2013 Lars Simon Zehnder
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

#' Finmix `mcmcoutputpermfixhierpost` class 
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
#' Note this class inherits all slots from its parent classes.  
#' 
#' @exportClass mcmcoutputpermfixhierpost
#' @describeIn mcmcoutputperm_class
#' @seealso 
#' * [mcmcoutputfixhierpost][mcmcoutput_class] for the parent class
#' * [mcmcpermfix][mcmcperm_class] for the parent class
#' * [mcmcpermute()] for performing permutation of MCMC samples
.mcmcoutputpermfixhierpost <- setClass("mcmcoutputpermfixhierpost",
  contains = c(
    "mcmcpermfixhier",
    "mcmcpermfixpost",
    "mcmcoutputfixhierpost"
  ),
  validity = function(object) {
    ## else: OK
    TRUE
  }
)

#' Initializer of the `mcmcoutputpermfixhier` class
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
#' @param logperm A named list containing the mixture log-likelihood, the 
#'   prior log-likelihood, and the complete data posterior log-likelihood 
#'   for the permuted MCMC samples.
#' @param hyperperm A named list containing the (permuted) parameters of the 
#'   hierarchical prior.
#' @param postperm A named list containing a named list `par` with array(s) of 
#'   parameters from the posterior density. 
#' @exportMethod initialize
#' @keywords internal
#' 
#' @seealso 
#' * [Classes_Details] for details of class definitions, and 
#' * [setOldClass] for the relation to S3 classes
setMethod(
  "initialize", "mcmcoutputpermfixhierpost",
  function(.Object, mcmcoutput, Mperm = integer(),
           parperm = list(), logperm = list(),
           hyperperm = list(), postperm = list()) {
    .Object@M <- mcmcoutput@M
    .Object@burnin <- mcmcoutput@burnin
    .Object@ranperm <- mcmcoutput@ranperm
    .Object@par <- mcmcoutput@par
    .Object@log <- mcmcoutput@log
    .Object@hyper <- mcmcoutput@hyper
    .Object@post <- mcmcoutput@post
    .Object@model <- mcmcoutput@model
    .Object@prior <- mcmcoutput@prior
    .Object@Mperm <- Mperm
    .Object@parperm <- parperm
    .Object@logperm <- logperm
    .Object@hyperperm <- hyperperm
    .Object@postperm <- postperm
    .Object
  }
)

#' Shows a summary of an `mcmcoutputpermfixhierpost` object.
#' 
#' @description
#' Calling [show()] on an `mcmcoutputpermfixhierpost` object gives an overview 
#' of the `mcmcoutputpermfixhierpost` object.
#' 
#' @param object An `mcmcoutputpermfixhierpost` object.
#' @returns A console output listing the slots and summary information about
#'   each of them. 
#' @exportMethod show
#' @describeIn mcmcoutputpermfixhierpost-class
setMethod(
  "show", "mcmcoutputpermfixhierpost",
  function(object) {
    cat("Object 'mcmcoutputperm'\n")
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
      "     post        : List of",
      length(object@post), "\n"
    )
    cat("     Mperm       :", object@Mperm, "\n")
    cat(
      "     parperm     : List of",
      length(object@parperm), "\n"
    )
    cat(
      "     logperm     : List of",
      length(object@logperm), "\n"
    )
    cat(
      "     postperm    : List of",
      length(object@postperm), "\n"
    )
    cat(
      "     hyperperm   : List of",
      length(object@hyperperm), "\n"
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
#' parameters are plotted together with `K-1` weights. If a hierarchical prior 
#' has been used in sampling its parameters are plotted alongside the other 
#' parameters.
#' 
#' @param x An `mcmcoutputpermfixhierpost` object containing all sampled values.
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
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2, 
#'                  indicfix = TRUE)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc()
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
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
    x = "mcmcoutputpermfixhierpost",
    dev = "ANY",
    lik = "ANY",
    col = "ANY"
  ),
  function(x, dev = TRUE, lik = 1, col = FALSE, ...) {
    dist <- x@model@dist
    if (lik %in% c(0, 1)) {
      if (dist == "poisson") {
        .permtraces.Poisson.Hier(x, dev)
      } else if (dist == "binomial") {
        .permtraces.Binomial(x, dev)
      } else if (dist == "exponential") {
        .permtraces.Exponential(x, dev)
      } else if (dist == "normal") {
        .permtraces.Normal.Hier(x, dev)
      } else if (dist == "student") {
        .permtraces.Student.Hier(x, dev)
      } else if (dist == "normult") {
        .permtraces.Normult.Hier(x, dev, col)
      } else if (dist == "studmult") {
        .permtraces.Studmult.Hier(x, dev, col)
      }
    }
    if (lik %in% c(1, 2)) {
      ## log ##
      .permtraces.Log(x, dev)
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
#' Note, this method is so far only implemented for Poisson and Binomial 
#' mixture distributions. 
#' 
#' @param x An `mcmcoutputpermfixhierpost` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Histograms of the MCMC samples.
#' @exportMethod plotHist
#' 
#' @examples 
#' \dontrun{
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2, 
#'                  indicfix = TRUE)
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
    x = "mcmcoutputpermfixhierpost",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    dist <- x@model@dist
    if (dist == "poisson") {
      .permhist.Poisson.Hier(x, dev)
    } else if (dist == "binomial") {
      .permhist.Binomial(x, dev)
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
#' Note, this method is so far only implemented for mixtures of Poisson or 
#' Binomial distributions.
#' 
#' @param x An `mcmcoutputpermfixhierpost` object containing all sampled values.
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
#' f_mcmc <- mcmc()
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' f_outputperm <- mcmcpermute(f_output)
#' plotDens(f_outputperm)
#' }
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
    x = "mcmcoutputpermfixhierpost",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    dist <- x@model@dist
    if (dist == "poisson") {
      .permdens.Poisson.Hier(x, dev)
    } else if (dist == "binomial") {
      .permdens.Binomial(x, dev)
    }
  }
)

#' Plot point processes of the component parameters
#' 
#' @description 
#' Calling [plotPointProc()] plots point processes of the sampled component 
#' parameters from MCMC sampling.  
#' 
#' Note, this method is so far only implemented for mixtures of Poisson or 
#' Binomial distributions.
#' 
#' @param x An `mcmcoutputpermfixhierpost` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Point process of the MCMC samples.
#' @exportMethod plotPointProc
#' 
#' @examples 
#' \dontrun{
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2, 
#'                  indicfix = TRUE)
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
    x = "mcmcoutputpermfixhierpost",
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
#' Note, this method is so far only implemented for mixtures of Poisson or 
#' Binomial distributions. 
#' 
#' @param x An `mcmcoutputpermfixhierpost` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Sampliing representations of the MCMC samples.
#' @exportMethod plotSampRep
#' 
#' @examples 
#' \dontrun{
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2, 
#'                  indicfix = TRUE)
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
    x = "mcmcoutputpermfixhierpost",
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
#' Note, this method is so far only implemented for Poisson and Binomial 
#' mixture distributions.
#' 
#' @param x An `mcmcoutputpermfixhierpost` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Posterior densities of the MCMC samples.
#' @exportMethod plotPostDens
#' 
#' @examples 
#' \dontrun{
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2, 
#'                  indicfix = TRUE)
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
    x = "mcmcoutputpermfixhierpost",
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