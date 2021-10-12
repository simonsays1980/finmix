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

#' Finmix `mcmcoutputfixhierpost` class 
#' 
#' @description 
#' This class inherits from the `mcmcoutputfixhier` class and adds posterior 
#' density parameters to the MCMC sampling output. The storage of posterior 
#' parameters is controlled by the slot `storepost` in the [mcmc-class] 
#' class. If set to `TRUE` posterior parameters are stored in the output of the 
#' MCMC sampling.
#' 
#' @slot post A named list containing a named list `par` with arrays for the 
#'   posterior density parameters.
#' @exportClass mcmcoutputfixhierpost
#' @rdname mcmcoutputfixhierpost-class
#' 
#' @seealso 
#' * [mcmcoutputfixhier-class] for the parent class
#' * [mixturemcmc()] for performing MCMC sampling
#' * [mcmc-class] for the class defining the MCMC hyper-parameters 
#' * [mcmc()] for the `mcmc` class constructor
.mcmcoutputfixhierpost <- setClass("mcmcoutputfixhierpost",
  representation(post = "list"),
  contains = c("mcmcoutputfixhier"),
  validity = function(object) {
    ## else: OK
    TRUE
  },
  prototype(post = list())
)

#' Shows a summary of an `mcmcoutputfixhierpost` object.
#' 
#' @description
#' Calling [show()] on an `mcmcoutputfixhierpost` object gives an overview 
#' of the `mcmcoutputfixhierpost` object.
#' 
#' @param object An `mcmcoutputfixhierpost` object.
#' @returns A console output listing the slots and summary information about
#'   each of them. 
#' @exportMethod show
#' @noRd
setMethod(
  "show", "mcmcoutputfixhierpost",
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
      "     post        : List of",
      length(object@post), "\n"
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
#' set to `1`. If `lik` is set to `0` the parameters of the components and the 
#' posterior parameters are plotted together with `K-1` weights.
#' 
#' Note that this method calls the equivalent method from the parent class 
#' `mcmcoutputfixhier`. 
#' 
#' @param x An `mcmcoutputfixhierpost` object containing all sampled values.
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
    x = "mcmcoutputfixhierpost",
    dev = "ANY",
    lik = "ANY",
    col = "ANY"
  ),
  function(x, dev = TRUE, lik = 1, col = FALSE, ...) {
    ## Call method 'plot()' from 'mcmcoutputfixhier'
    callNextMethod(x, dev, lik, col, ...)
  }
)

#' Plot histograms of the parameters and weights
#' 
#' @description 
#' Calling [plotHist()] plots histograms of the sampled parameters and weights 
#' from MCMC sampling.More specifically, all component parameters, `K-1` of the 
#' weights and the posterior parameters are considered in the histogram plots. 
#' 
#' Note that this method calls the equivalent method from the parent class 
#' `mcmcoutputfixhier`.
#' 
#' @param x An `mcmcoutputfixhierpost` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Histograms of the MCMC samples.
#' @exportMethod plotHist
#' @noRd
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
    x = "mcmcoutputfixhierpost",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    ## Call 'plotHist()' from 'mcmcoutputfixhier'
    callNextMethod(x, dev, ...)
  }
)

#' Plot densities of the parameters and weights
#' 
#' @description 
#' Calling [plotDens()] plots densities of the sampled parameters and weights 
#' from MCMC sampling.More specifically, all component parameters, `K-1` of the 
#' weights and the posterior parameters are considered in the density plots. 
#' 
#' Note that this method calls the equivalent method from the parent class 
#' `mcmcoutputfixhier`. 
#' 
#' @param x An `mcmcoutputfixhierpost` object containing all sampled values.
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
    x = "mcmcoutputfixhierpost",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    ## Call 'plotDens()' from 'mcmcoutputfixhier'
    callNextMethod(x, dev, ...)
  }
)

#' Plot point processes of the component parameters
#' 
#' @description 
#' Calling [plotPointProc()] plots point processes of the sampled component 
#' parameters from MCMC sampling.  
#' 
#' Note, this methid calls the equivalent method from the parent class 
#' `mcmcoutputfixhier`.
#' 
#' @param x An `mcmcoutputfixhierpost` object containing all sampled values.
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
    x = "mcmcoutputfixhierpost",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    ## Call 'plotPointProc()' from 'mcmcoutputfixhier'
    callNextMethod(x, dev, ...)
  }
)

#' Plot sampling representations for the component parameters.
#' 
#' @description 
#' Calling [plotSampRep()] plots sampling representations of the sampled 
#' component parameters from MCMC sampling.  
#' 
#' Note, this method calls the equivalent method of the parent class 
#' `mcmcoutputfixhier`.
#' 
#' @param x An `mcmcoutputfixhierpost` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Sampling representation of the MCMC samples.
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
    x = "mcmcoutputfixhierpost",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    ## Call 'plotSampRep()' from 'mcmcoutputfixhier'
    callNextMethod(x, dev, ...)
  }
)

#' Plot posterior densities of the component parameters
#' 
#' @description 
#' Calling [plotPostDens()] plots posterior densities of the sampled component 
#' parameters from MCMC sampling, if the number of components is two. 
#' 
#' Note, this method calls the equivalent method of the parent class 
#' `mcmcoutputfixhier`.
#' 
#' @param x An `mcmcoutputfixhierpost` object containing all sampled values.
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
    x = "mcmcoutputfixhierpost",
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
#' Note that this method calls the equivalent method from the parent class and 
#' adds the sub-chains for the posterior density parameters.
#' 
#' @param object An `mcmcoutput` object containing all sampled values.
#' @param index An array specifying the extraction of the sub-chain.
#' @return An `mcmcoutput` object containing the values from the sub-chain.
#' @noRd
#' @export subseq
setMethod(
  "subseq", signature(
    object = "mcmcoutputfixhierpost",
    index = "array"
  ),
  function(object, index) {
    ## TODO: Check arguments via .validObject ##
    dist <- object@model@dist
    ## Call 'subseq()' from 'mcmcoutputfixhier'
    object <- callNextMethod(object, index)
    ## post ##
    if (dist == "poisson") {
      .subseq.Poisson.Post(object, index)
    } else if (dist == "binomial") {
      .subseq.Binomial.Mcmcoutputfixpost(object, index)
    } else if (dist %in% c("normal", "student")) {
      .subseq.Norstud.Mcmcoutputfixpost(object, index)
    } else if (dist %in% c("normult", "studmult")) {
      .subseq.Normultstud.Mcmcoutputfixpost(object, index)
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
    object = "mcmcoutputfixhierpost",
    index = "array"
  ),
  function(object, index) {
    if (object@model@K == 1) {
      return(object)
    } else {
      ## Check arguments, TODO: .validObject ##
      dist <- object@model@dist
      ## Call 'swapElements()' from 'mcmcoutputfixhier'
      object <- callNextMethod(object, index)
      if (dist == "poisson") {
        .swapElements.Poisson.Post(object, index)
      } else if (dist == "binomial") {
        .swapElements.Binomial.Mcmcoutputfixpost(object, index)
      } else if (dist %in% c("normal", "student")) {
        .swapElements.Norstud.Mcmcoutputfixpost(object, index)
      } else if (dist %in% c("normult", "studmult")) {
        .swapElements.Normultstud.Mcmcoutputfixpost(object, index)
      }
    }
  }
)
