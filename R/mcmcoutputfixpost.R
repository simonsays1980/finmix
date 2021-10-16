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

#' Finmix `mcmcoutput` class for fixed indicators and posterior parameters
#' 
#' @description 
#' The `mcmcoutputfixpost` class inherits from the `mcmcoutputfix` class and 
#' adds a slot to store the parameters of the posterior distribution from which 
#' the component parameters are drawn. The storage of posterior parameters is 
#' controlled by the slot `storepost` in the [mcmc-class] class. If set 
#' to `TRUE` posterior parameters are stored in the output of the MCMC sampling. 
#' 
#' @slot post A named list containing a list `par` that contains the posterior 
#'   parameters as named arrays. 
#' @exportClass mcmcoutputfixpost
#' @rdname mcmcoutputfixpost-class
#' @keywords internal
#' 
#' @seealso 
#' * [mcmcoutputfix-class] for the parent class
#' * [mcmcoutputpost-class] for the corresponding class for unknown 
#'   indicators.
#' * [mcmc-class] for the class defining the MCMC hyper-parameters
#' * [mcmc()] for the constructor of the [mcmc-class] class
.mcmcoutputfixpost <- setClass("mcmcoutputfixpost",
  representation(post = "list"),
  contains = c("mcmcoutputfix"),
  validity = function(object) {
    ## else: OK
    TRUE
  },
  prototype(post = list())
)

#' Shows a summary of an `mcmcoutputfixpost` object.
#' 
#' @description
#' Calling [show()] on an `mcmcoutputfixpost` object gives an overview 
#' of the `mcmcoutputfixpost` object.
#' 
#' @param object An `mcmcoutputfixpost` object.
#' @returns A console output listing the slots and summary information about
#'   each of them. 
#' @exportMethod show
#' @noRd
setMethod(
  "show", "mcmcoutputfixpost",
  function(object) {
    cat("Object 'mcmcoutputfixpost\n")
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
#' `mcmcoutputfix`. 
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
#' f_mcmc <- mcmc()
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
    x = "mcmcoutputfixpost",
    dev = "ANY",
    lik = "ANY",
    col = "ANY"
  ),
  function(x, dev = TRUE, lik = 1, col = FALSE, ...) {
    ## Call 'plot()' from 'mcmcoutputfix
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
#' `mcmcoutputfix`.
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
#' f_mcmc <- mcmc()
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
    x = "mcmcoutputfixpost",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    ## Call 'plotHist()' from 'mcmcoutputfix'
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
#' Note that this methid calls the equivalent method from the parent class 
#' `mcmcoutputfix`. 
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
#' f_mcmc <- mcmc()
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
    x = "mcmcoutputfixpost",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    ## Call 'plotDens()' from 'mcmcoutputfix'
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
#' `mcmcoutputfix`.
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
#' f_mcmc <- mcmc()
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
    x = "mcmcoutputfixpost",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    ## Call 'plotPointProc()' from 'mcmcoutputfix'
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
#' `mcmcoutputfix`.
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
#' f_mcmc <- mcmc()
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
    x = "mcmcoutputfixpost",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    ## Call 'plotSampRep()' from 'mcmcoutputfix'
    callNextMethod(x, dev, ...)
  }
)

#' Plot posterior densities of the component parameters
#' 
#' @description 
#' Calling [plotPostDens()] plots posterior densities of the sampled component 
#' parameters from MCMC sampling, if the number of components is two. 
#' 
#' Note, this methid calls the equivalent method of the parent class 
#' `mcmcoutputfix`.
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
#' f_mcmc <- mcmc()
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
    x = "mcmcoutputfixpost",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    ## Call 'plotPostDens()' from 'mcmcoutputfix'
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
#' @exportMethod subseq
#' @noRd
setMethod(
  "subseq", signature(
    object = "mcmcoutputfixpost",
    index = "array"
  ),
  function(object, index) {
    ## Call 'subseq()' from 'mcmcoutputfix'
    callNextMethod(object, index)
    dist <- object@model@dist
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
#' @exportMethod swapElements
#' @noRd 
setMethod(
  "swapElements", signature(
    object = "mcmcoutputfixpost",
    index = "array"
  ),
  function(object, index) {
    if (object@model@K == 1) {
      return(object)
    } else {
      ## Call method 'swapiElements()' from 'mcmcoutputfix'
      object <- callNextMethod()
      dist <- object@model@dist
      if (dist == "poisson") {
        .swapElements.Poisson(object, index)
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

#' Getter method of `mcmcoutputfixpost` class.
#' 
#' Returns the `post` slot.
#' 
#' @param object An `mcmcoutputfixpost` object.
#' @returns The `post` slot of the `object`.
#' @exportMethod getPost
#' @noRd
#' 
#' @examples 
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
#' # Get the slot.
#' getPost(f_output)
#' 
#' @seealso 
#' * [mcmcoutput-class] for the class definition
#' * [mixturemcmc()] for performing MCMC sampling
setMethod(
  "getPost", "mcmcoutputfixpost",
  function(object) {
    return(object@post)
  }
)

## No setters as users are not intended to manipulate ##
## this object. ##

#' Generates sub-chains from Poisson MCMC samples
#' 
#' @description 
#' For internal usage only. This function generates sub-chains from an 
#' `mcmcoutput` object with a Poisson `model` by defining an `index` array 
#' specifying how extraction of sub-samples should be performed. Has errors for 
#' some `mcmcoutput` sub-classes. Note, this method is only supplementing the 
#' method from the parent class by adding the sub-chains for the parameters of 
#' the posterior density. 
#' 
#' 
#' @param obj An `mcmcoutputfixpost` object containing all MCMC samples.
#' @param index An array specifying the extraction of sub-samples.
#' @return An `mcmcoutputfixpost` object containing sub-chains.
#' @noRd
#' 
#' @seealso 
#' * [subseq()] for the calling method
".subseq.Poisson.Post" <- function(obj, index) {
  if (obj@model@K == 1) {
    obj@post$par$a <- array(obj@post$par$a[index],
      dim = c(obj@M, 1)
    )
    obj@post$par$b <- array(obj@post$par$b[index],
      dim = c(obj@M, 1)
    )
  } else {
    obj@post$par$a <- obj@post$par$a[index, ]
    obj@post$par$b <- obj@post$par$b[index, ]
  }
  return(obj)
}

#' Swaps elements in Poisson MCMC samples.
#' 
#' @description 
#' For internal usage only. This function swaps elements for each row in an 
#' `mcmcoutput` object with Poisson MCMC samples. Calls the C++-function 
#' `swap_cc()`. Note that this function only complements the equivalent 
#' functionality of the parent class by also swapping the posterior density 
#' parameters. 
#' 
#' @param obj An `mcmcoutput` object containing all MCMC samples.
#' @param index An array specifying the element swapping.
#' @return An `mcmcoutput` object with swapped elements.
#' @noRd
#' 
#' @seealso 
#' * [swapElements()] for the calling method
".swapElements.Poisson.Post" <- function(obj, index) {
  ## Rcpp::export 'swap_cc'
  obj@post$par$a <- swap_cc(obj@post$par$a, index)
  obj@post$par$b <- swap_cc(obj@post$par$b, index)
  return(obj)
}

#' Generates sub-chains from Binomial MCMC samples
#' 
#' @description 
#' For internal usage only. This function generates sub-chains from an 
#' `mcmcoutput` object with a Binomial `model` by defining an `index` array 
#' specifying how extraction of sub-samples should be performed. Has errors for 
#' some `mcmcoutput` sub-classes. Note, this method is only supplementing the 
#' method from the parent class by adding the sub-chains for the parameters of 
#' the posterior density. 
#' 
#' 
#' @param obj An `mcmcoutputfixpost` object containing all MCMC samples.
#' @param index An array specifying the extraction of sub-samples.
#' @return An `mcmcoutputfixpost` object containing sub-chains.
#' @noRd
#' 
#' @seealso 
#' * [subseq()] for the calling method
".subseq.Binomial.Mcmcoutputfixpost" <- function(obj, index) {
  if (obj@model@K == 1) {
    obj@post$par$a <- array(obj@post$par$a[index],
      dim = c(obj@M, 1)
    )
    obj@post$par$b <- array(obj@post$par$b[index],
      dim = c(obj@M, 1)
    )
  } else {
    obj@post$par$a <- obj@post$par$a[index, ]
    obj@post$par$b <- obj@post$par$b[index, ]
  }
  return(obj)
}

#' Swaps elements in Binomial MCMC samples.
#' 
#' @description 
#' For internal usage only. This function swaps elements for each row in an 
#' `mcmcoutput` object with Binomial MCMC samples. Calls the C++-function 
#' `swap_cc()`. Note that this function only complements the equivalent 
#' functionality of the parent class by also swapping the posterior density 
#' parameters. 
#' 
#' @param obj An `mcmcoutput` object containing all MCMC samples.
#' @param index An array specifying the element swapping.
#' @return An `mcmcoutput` object with swapped elements.
#' @noRd
#' 
#' @seealso 
#' * [swapElements()] for the calling method
".swapElements.Binomial.Mcmcoutputfixpost" <- function(obj, index) {
  ## Rcpp::export 'swap_cc'
  obj@post$par$a <- swap_cc(obj@post$par$a, index)
  obj@post$par$b <- swap_cc(obj@post$par$b, index)
  return(obj)
}

#' Generates sub-chains from Normal or Student-t MCMC samples
#' 
#' @description 
#' For internal usage only. This function generates sub-chains from an 
#' `mcmcoutput` object with a Normal or Student-t `model` by defining an 
#' `index` array specifying how extraction of sub-samples should be performed. 
#' Has errors for some `mcmcoutput` sub-classes. Note, this method is only 
#' supplementing the method from the parent class by adding the sub-chains for 
#' the parameters of the posterior density. 
#' 
#' 
#' @param obj An `mcmcoutputfixpost` object containing all MCMC samples.
#' @param index An array specifying the extraction of sub-samples.
#' @return An `mcmcoutputfixpost` object containing sub-chains.
#' @noRd
#' 
#' @seealso 
#' * [subseq()] for the calling method
".subseq.Norstud.Mcmcoutputfixpost" <- function(obj, index) {
  if (obj@model@K == 1) {
    obj@post$par$mu$b <- array(obj@post$par$mu$b[index],
      dim = c(obj@M, 1)
    )
    obj@post$par$mu$B <- array(obj@post$par$mu$B[index],
      dim = c(obj@M, 1)
    )
    obj@post$par$sigma$c <- array(obj@post$par$sigma$c[index],
      dim = c(obj@M, 1)
    )
    obj@post$par$sigma$C <- array(obj@post$par$sigma$C[index],
      dim = c(obj@M, 1)
    )
  } else {
    obj@post$par$mu$b <- obj@post$par$mu$b[index, ]
    obj@post$par$mu$b <- obj@post$par$mu$B[index, ]
    obj@post$par$sigma$c <- obj@post$par$sigma$c[index, ]
    obj@post$par$sigma$C <- obj@post$par$sigma$C[index, ]
  }
  return(obj)
}

#' Swaps elements in Normal or Student-t MCMC samples.
#' 
#' @description 
#' For internal usage only. This function swaps elements for each row in an 
#' `mcmcoutput` object with Normal or Student-t MCMC samples. Calls the 
#' C++-function `swap_cc()`. Note that this function only complements the 
#' equivalent functionality of the parent class by also swapping the posterior 
#' density parameters. 
#' 
#' @param obj An `mcmcoutput` object containing all MCMC samples.
#' @param index An array specifying the element swapping.
#' @return An `mcmcoutput` object with swapped elements.
#' @noRd
#' 
#' @seealso 
#' * [swapElements()] for the calling method
".swapElements.Norstud.Mcmcoutputfixpost" <- function(obj, index) {
  ## Rcpp::export 'swap_cc'
  obj@post$par$mu$b <- swap_cc(obj@post$par$mu$b, index)
  obj@post$par$mu$B <- swap_cc(obj@post$par$mu$B, index)
  obj@post$par$sigma$c <- swap_cc(obj@post$par$sigma$c, index)
  obj@post$par$sigma$C <- swap_cc(obj@post$par$sigma$C, index)
  return(obj)
}

#' Generates sub-chains from multivariate Normal or Student-t MCMC samples
#' 
#' @description 
#' For internal usage only. This function generates sub-chains from an 
#' `mcmcoutput` object with a multivariate Normal or Student-t `model` by 
#' defining an `index` array specifying how extraction of sub-samples should be 
#' performed. Has errors for some `mcmcoutput` sub-classes. Note, this method 
#' is only supplementing the method from the parent class by adding the 
#' sub-chains for the parameters of the posterior density. 
#' 
#' 
#' @param obj An `mcmcoutputfixpost` object containing all MCMC samples.
#' @param index An array specifying the extraction of sub-samples.
#' @return An `mcmcoutputfixpost` object containing sub-chains.
#' @noRd
#' 
#' @seealso 
#' * [subseq()] for the calling method
".subseq.Normultstud.Mcmcoutputfixpost" <- function(obj, index) {
  if (obj@model@K == 1) {
    obj@post$par$mu$b <- obj@post$par$mu$b[index, ]
    obj@post$par$mu$B <- obj@post$par$mu$B[index, ]
    obj@post$par$sigma$c <- obj@post$par$sigma$c[index, ]
    obj@post$par$sigma$C <- obj@post$par$sigma$C[index, ]
  } else {
    obj@post$par$mu$b <- obj@post$par$mu$b[index, , ]
    obj@post$par$mu$B <- obj@post$par$mu$B[index, , ]
    obj@post$par$sigma$c <- obj@post$par$sigma$c[index, , ]
    obj@post$par$sigma$C <- obj@post$par$sigma$C[index, , ]
  }
  return(obj)
}

#' Swaps elements in multivariate Normal or Student-t MCMC samples.
#' 
#' @description 
#' For internal usage only. This function swaps elements for each row in an 
#' `mcmcoutput` object with multivariate Normal or Student-t MCMC samples. 
#' Calls the C++-function `swap_cc()`. Note that this function only complements 
#' the equivalent functionality of the parent class by also swapping the 
#' posterior density parameters. 
#' 
#' @param obj An `mcmcoutput` object containing all MCMC samples.
#' @param index An array specifying the element swapping.
#' @return An `mcmcoutput` object with swapped elements.
#' @noRd
#' 
#' @seealso 
#' * [swapElements()] for the calling method
".swapElements.Normultstud.Mcmcoutputfixpost" <- function(obj, index) {
  ## Rcpp::export 'swap_3d_cc'
  obj@post$par$mu$b <- swap_cc(obj@post$par$mu$b, index)
  obj@post$par$mu$B <- swap_3d_cc(obj@post$par$mu$B, index)
  obj@post$par$sigma$c <- swap_cc(obj@post$par$sigma$c, index)
  obj@post$par$sigma$C <- swap_3d_cc(obj@post$par$sigma$C, index)
  return(obj)
}
