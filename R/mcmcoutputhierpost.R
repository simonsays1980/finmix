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

#' Finmix `mcmcoutputhierpost` class
#' 
#' @description 
#' This class stores samples from bayesian estimation with hierarchical prior 
#' and unknown indicators. It inherits from `mcmcoutputhier` and adds to it a 
#' slot to store the parameters from the posterior density. For a model with 
#' unknown indicators the slot `@@indicfix` in the `model` object specifying 
#' the finite mixture model must be set to `FALSE` (default). Sampling with a 
#' hierarchical prior is activated by setting the slot `@@hier` in the `prior` 
#' object to `TRUE` (default). Finally, to store parameters for the posterior 
#' density the hyper-parameter `storepost` in the `mcmc` object must be set to 
#' `TRUE` (default).
#' 
#' @slot post A named list containing a named list `par` that contains arrays 
#'   storing the sampled posterior density parameters.
#' @exportClass mcmcoutputhierpost
#' @rdname mcmcoutputhierpost-class
#' 
#' @seealso 
#' * [mcmcoutputhier-class] for the parent class
#' * [prior-class] for the class specifying the prior distribution 
#' * [prior()] for the `prior` class constructor
#' * [priordefine()] for the advanced `prior` class constructor 
#' * [mcmc-class] for the class defining the hyper-parameters
#' * [mcmc()] for the `mcmc` class constructor
#' * [model-class] for the `model` class definition
#' * [model()] for the `model` class constructor
.mcmcoutputhierpost <- setClass("mcmcoutputhierpost",
  representation(post = "list"),
  contains = c("mcmcoutputhier"),
  validity = function(object) {
    ## else: OK
    TRUE
  },
  prototype(post = list())
)

#' Finmix `mcmcoutput` class
#' 
#' @description 
#' The `mcmcoutput` class stores all MCMC samples and corresponding information.
#' 
#' @details  
#' Calling [mixturemcmc()] on appropriate input arguments performs MCMC 
#' sampling and returns an `mcmcoutput` object that stores all samples and 
#' corresponding information like hyper-parameters, the finite mixture model 
#' specified in a `model` object and the `prior` that specifies the prior 
#' distribution. All slots are listed below. Note that not all slots must be 
#' available in a object of class `mcmcoutput`. Some slots get only occupied, 
#' if a hierarchical prior had been used in MCMC sampling, or if posterior 
#' samples should be stored. Furthermore, the slots also look different, if 
#' MCMC sampling had been performed for a model with fixed indicators (see for 
#' subclasses for example [mcmcoutputfix-class], [mcmcoutputbase-class],
#' [mcmcoutputhier-class] or [mcmcoutputpost-class]). 
#' 
#' The class `mcmcoutput` is a class union and includes all classes that 
#' define objects to store MCMC samples and is used to dispatch methods for 
#' `mcmcoutput` objects. For the user this detail is not important, 
#' especially as this class has no exported constructor. Objects are solely 
#' constructed internally within the function [mixturemcmc()]. 
#' 
#' ## Class methods 
#' 
#' This class comes along with a couple of methods that should give the user 
#' some comfort in handling the MCMC sampling results. There are no setters for 
#' this class as the slots are only set internally. 
#' 
#' ### Show 
#' * `show()` shows a short summary of the object's slots.
#' 
#' ### Getters
#' * `getM()` returns the `M` slot.
#' * `getBurnin()` returns the `burnin` slot.
#' * `getRanperm()` returns the `ranperm` slot.
#' * `getPar()` returns the `par` slot.
#' * `gteWeight()` returns the `weight` slot, if available.
#' * `getLog()` returns the `log` slot.
#' * `getEntropy()` returns the `entropy` slot, if available.
#' * `getHyper()` returns the `hyper` slot, if available.
#' * `getPost()` returns the `post` slot, if available.
#' * `getST()` returns the `ST` slot, if available.
#' * `getS()` returns the `S` slot, if available.
#' * `getNK()` returns the `NK` slot, if available.
#' * `getClust()` returns the `clust` slot, if available.
#' * `getModel()` returns the `model` slot. 
#' * `getPrior()` returns the `prior` slot.   
#' 
#' ### Plotting
#' Plotting functionality for the `mcmcoutput` helps the user to inspect MCMC 
#' results.
#' 
#' * `plotTraces()` plots traces of MCMC samples. See [plotTraces()] for 
#'   further information.
#' * `plotHist()` plots histograms of parameters and weights. See [plotHist()] 
#'   for further information.
#' * `plotDens()` plots densities of parameters and weights. See [plotDens()] 
#'   for further information.
#' * `plotPointProc()` plots the point process of component parameters. See 
#'   [plotPointProc] for further information.
#' * `plotSampRep()` plots the sampling representation of component parameters. 
#'   See [plotSampRep()] for further information. 
#' * `plotPostDens()` plots the posterior density of component parameters. Note 
#'   that this function can only be applied for mixtures of two components. See 
#'   [plotPostDens()] for further information.  
#'
#' ## Slots
#' * `M` An integer defining the number of iterations in MCMC sampling. 
#' * `burnin` An integer defining the number of iterations in the burn-in 
#'   phase of MCMC sampling. These number of sampling steps are not stored 
#'   in the output.
#' * `ranperm` A logical indicating, if MCMC sampling has been performed 
#'   with random permutations of components.
#' * `par` A named list containing the sampled component parameters. 
#' * `weight` An `array` of dimension `M x K` containing the sampled 
#'   weights.
#' * `log` A named list containing the values of the mixture log-likelihood, 
#'   mixture prior log-likelihood, and the complete data posterior 
#'   log-likelihood.
#' * `hyper` A list storing the sampled parameters from the hierarchical 
#'   prior. 
#' * `post` A named list containing a list `par` that contains the posterior 
#'   parameters as named arrays. 
#' * `entropy` An `array` of dimension `M x 1` containing the entropy 
#'   for each MCMC draw.
#' * `ST` An `array` of dimension `M x 1` containing all MCMC states, 
#'   for the last observation in slot `y` of the `fdata` object passed in to 
#'   [mixturemcmc()] where a state is defined for non-Markov models as the 
#'   last indicator of this observation. 
#' * `S` An `array` of dimension `N x storeS` containing the last 
#'   `storeS` indicators sampled. `storeS` is defined in the slot `storeS` of 
#'   the `mcmc` object passed into [mixturemcmc()].
#' * `NK` An `array` of dimension `M x K` containing the number of 
#'   observations assigned to each component for each MCMC draw.
#' * `clust` An `array` of dimension `N x 1` containing the recent 
#'   indicators defining the last "clustering" of observations into the 
#'   mixture components.
#' * `model` The `model` object that specifies the finite mixture model for 
#'   which MCMC sampling has been performed. 
#' * `prior` The `prior` object defining the prior distributions for the 
#'   component parameters that has been used in MCMC sampling.
#'   
#' @exportClass mcmcoutput
#' @rdname mcmcoutput-class
#' 
#' @seealso 
#' * [mcmcoutputperm-class] for the corresponding class defined for relabeled 
#'   MCMC samples
#' * [mcmcoutputfix-class] for the `mcmcoutput` sub-class for models with 
#'   fixed indicators
#' * [mcmcoutputbase-class] for the `mcmcoutput` sub-class for models with 
#'   unknown indicators
#' * [mcmcoutputhier-class] for the `mcmcoutput` sub-class for MCMC samples 
#'   with hierarchical priors
#' * [mcmcoutputpost-class] for the `mcmcoutput` sub-class for MCMC samples 
#'   with stored posterior density parameters
setClassUnion(
  "mcmcoutput",
  c(
    "mcmcoutputfix",
    "mcmcoutputfixhier",
    "mcmcoutputfixpost",
    "mcmcoutputfixhierpost",
    "mcmcoutputbase",
    "mcmcoutputhier",
    "mcmcoutputpost",
    "mcmcoutputhierpost"
  )
)

#' Shows a summary of an `mcmcoutputhierpost` object.
#' 
#' @description
#' Calling [show()] on an `mcmcoutputhierpost` object gives an overview 
#' of the `mcmcoutputhierpost` object.
#' 
#' @param object An `mcmcoutputhierpost` object.
#' @returns A console output listing the slots and summary information about
#'   each of them. 
#' @exportMethod show
#' @keywords internal
setMethod(
  "show", "mcmcoutputhierpost",
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
#' set to `1`.s If `lik` is set to `0` the parameters of the components and the 
#' posterior parameters are plotted together with `K-1` weights. 
#' 
#' Note that this method calls the equivalent method from the parent class.
#' 
#' @param x An `mcmcoutputhierpost` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param lik An integer indicating, if the log-likelihood traces should be 
#'   plotted (default). If set to `0` the traces for the parameters 
#'   and weights are plotted instead. 
#' @param col A logical indicating, if the plot should be colored.
#' @param ... Further arguments to be passed to the plotting function.
#' @return A plot of the traces of the MCMC samples.
#' @exportMethod plotTraces
#' @keywords internal
#' 
#' @examples 
#' \dontrun{
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
    x = "mcmcoutputhierpost",
    dev = "ANY",
    lik = "ANY",
    col = "ANY"
  ),
  function(x, dev = TRUE, lik = 1, col = FALSE, ...) {
    ## Call method 'plot()' from 'mcmcoutputhier'
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
#' Note, this method calls the equivalent method of the parent class. 
#' 
#' @param x An `mcmcoutputhierpost` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Histograms of the MCMC samples.
#' @exportMethod plotHist
#' @keywords internal
#' 
#' @examples 
#' \dontrun{
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2)
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
    x = "mcmcoutputhierpost",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    ## Call method 'plotHist()' from 'mcmcoutputhier'
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
#' Note that this method calls the equivalent method of the parent class.
#' 
#' @param x An `mcmcoutputhierpost` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Densities of the MCMC samples.
#' @exportMethod plotDens
#' @keywords internal
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
    x = "mcmcoutputhierpost",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    ## Call 'plotHist()' from 'mcmcoutputhier'
    callNextMethod(x, dev, ...)
  }
)

#' Plot point processes of the component parameters
#' 
#' @description 
#' Calling [plotPointProc()] plots point processes of the sampled component 
#' parameters from MCMC sampling.
#' 
#' Note, this method calls the equivalent method of the parent class. 
#' 
#' @param x An `mcmcoutputhierpost` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Point process of the MCMC samples.
#' @exportMethod plotPointProc
#' @keywords internal
#' 
#' @examples 
#' \dontrun{
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
    x = "mcmcoutputhierpost",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    ## Call 'plotPointProc()' from 'mcmcoutputhier'
    callNextMethod(x, dev, ...)
  }
)

#' Plot sampling representations for the component parameters.
#' 
#' @description 
#' Calling [plotSampRep()] plots sampling representations of the sampled 
#' component parameters from MCMC sampling.  
#' 
#' Note, this method calls the equivalent method from the parent class.
#' 
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Sampling representation of the MCMC samples.
#' @exportMethod plotSampRep
#' @keywords internal
#' 
#' @examples 
#' \dontrun{
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
    x = "mcmcoutputhierpost",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    ## Call 'plotSampRep()' from 'mcmcoutputhier'
    callNextMethod(x, dev, ...)
  }
)

#' Plot posterior densities of the component parameters
#' 
#' @description 
#' Calling [plotPostDens()] plots posterior densities of the sampled component 
#' parameters from MCMC sampling, if the number of components is two. 
#' 
#' Note, this method calls the equivalent method of the parent class.
#' 
#' @param x An `mcmcoutputhierpost` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Posterior densities of the MCMC samples.
#' @exportMethod plotPostDens
#' @keywords internal
#' 
#' @examples 
#' \dontrun{
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
    x = "mcmcoutputhierpost",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    ## Call 'plotPostDens()' from 'mcmcoutputhier'
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
#' adds to it the sub-chains for the parameters of the posterior density by 
#' calling a function from the `mcmcoutputfixpost` class.
#' 
#' @param object An `mcmcoutputhierpost` object containing all sampled values.
#' @param index An array specifying the extraction of the sub-chain.
#' @return An `mcmcoutputhierpost` object containing the values from the 
#'   sub-chain.
#' @exportMethod subseq
#' @keywords internal
setMethod(
  "subseq", signature(
    object = "mcmcoutputhierpost",
    index = "array"
  ),
  function(object, index) {
    ## Call 'subseq()' method from 'mcmcoutputhier'
    as(object, "mcmcoutputhier") <- callNextMethod(object, index)
    # Change owned slots #
    dist <- object@model@dist
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
#' @param object An `mcmcoutputhierpost` object containing the sampled values.
#' @param index An array specifying the extraction of the values.
#' @return An `mcmcoutputhierpost` object with swapped elements.
#' @exportMethod swapElements
#' @keywords internal 
setMethod(
  "swapElements", signature(
    object = "mcmcoutputhierpost",
    index = "array"
  ),
  function(object, index) {
    ## Check arguments, TODO: .validObject ##
    if (object@model@K == 1) {
      return(object)
    } else {
      dist <- object@model@dist
      ## Call method 'swapElements()' from 'mcmcoutputhier'
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

#' Getter method of `mcmcoutputhierpost` class.
#' 
#' Returns the `post` slot.
#' 
#' @param object An `mcmcoutputhierpost` object.
#' @returns The `post` slot of the `object`.
#' @exportMethod getPost
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
#' # Get the slot.
#' getPost(f_output)
#' 
#' @seealso 
#' * [mcmcoutput-class] for the class definition
#' * [mixturemcmc()] for performing MCMC sampling
setMethod(
  "getPost", "mcmcoutputhierpost",
  function(object) {
    return(object@post)
  }
)

## No setters as users are not intended to manipuate ##
## this object ##
