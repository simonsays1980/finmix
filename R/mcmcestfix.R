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

#' Finmix `mcmcestfix` class
#' 
#' @description 
#' This class stores the point estimators for component parameters and weights 
#' as well as corresponding information from MCMC sampling. Three point 
#' estimators are calculated: the maximum a posterior (MAP), the Bayesian 
#' maximum likelihood (BML) and the Identified ergodic average (IEAVG). See 
#' Fr\"uhwirth-Schnatter (2006) for detailed information about how these 
#' estimators are defined. 
#' 
#' @slot dist A character specifying the distribution family of the mixture 
#'   model used in MCMC sampling.
#' @slot K An integer specifying the number of components in the mixture model. 
#' @slot indicmod A character specifying the indicator model. At this moment 
#'   only a multinomial model can be chosen. 
#' @slot burnin An integer specifying the number of iterations in the burn-in 
#'   phase of MCMC sampling. 
#' @slot M An integer specifying the number of iterations to store in MCMC 
#'   sampling.
#' @slot ranperm A logical specifying, if random permutation has been used 
#'   during MCMC sampling. 
#' @slot relabel A character specifying the re-labeling algorithm used during 
#'   parameter estimation for the identified ergodic average. 
#' @slot map A named list containing the parameter estimates of the MAP. The 
#'   element `par` is a named list and contains the component parameters and 
#'   the element `weight` contains the weights. 
#' @slot bml A named list containing the parameter estimates of the BML. The 
#'   element `par` is a named list and contains the component parameters and 
#'   the element `weight` contains the weights. 
#' @slot ieavg A named list containing the parameter estimates of the IEAVG. The 
#'   element `par` is a named list and contains the component parameters and 
#'   the element `weight` contains the weights.
#' @slot sdpost A named list containing the standard deviations of the 
#'   parameter estimates from the posterior distributions.
#' @exportClass mcmcestfix
#' @rdname mcmcestfix-class
#' @keywords internal
#' 
#' @seealso
#' * [mcmcestind-class] for the equivalent class for models with 
#'   unknown indicators
#' * [mcmcestimate()] to calculate point estimates
.mcmcestfix <- setClass("mcmcestfix",
  representation(
    dist = "character",
    K = "integer",
    indicmod = "character",
    burnin = "integer",
    M = "integer",
    ranperm = "logical",
    relabel = "character",
    map = "list",
    bml = "list",
    ieavg = "list",
    sdpost = "list"
  ),
  validity = function(object) {
    ## else: OK
    TRUE
  },
  prototype(
    dist = character(),
    K = integer(),
    indicmod = character(),
    burnin = integer(),
    M = integer(),
    ranperm = logical(),
    relabel = character(),
    map = list(),
    bml = list(),
    ieavg = list(),
    sdpost = list()
  )
)

#' Shows a summary of an `mcmcestfix` object.
#' 
#' Calling [show()] on an `mcmcestfix` object gives an overview 
#' of the `mcmcestfix` object.
#' 
#' @param object An `mcmcestfix` object.
#' @returns A console output listing the slots and summary information about
#'   each of them. 
#' @exportMethod show
#' @keywords internal
setMethod(
  "show", "mcmcestfix",
  function(object) {
    cat("Object 'mcmcest'\n")
    cat("     dist        :", object@dist, "\n")
    cat("     K           :", object@K, "\n")
    cat(
      "     indicmod    :", object@indicmod,
      "\n"
    )
    cat("     M           :", object@M, "\n")
    cat("     burnin      :", object@burnin, "\n")
    cat("     ranperm     :", object@ranperm, "\n")
    cat("     relabel     :", object@relabel, "\n")
    cat(
      "     map         : List of",
      length(object@map), "\n"
    )
    cat(
      "     bml         : List of",
      length(object@bml), "\n"
    )
    cat(
      "     ieavg       : List of",
      length(object@ieavg), "\n"
    )
    cat(
      "     sdpost      : List of",
      length(object@sdpost), "\n"
    )
  }
)

#' Shows an advanced summary of an `mcmcestfix` object.
#' 
#' Calling [Summary()] on an `mcmcestfix` object gives an advanced overview 
#' of the `mcmcestfix` object.
#' 
#' Note, this method is so far only implemented for mixtures of Poisson 
#' distributions.
#' 
#' @param object An `mcmcestfix` object.
#' @returns A console output listing the formatted slots and summary 
#'   information about each of them. 
#' @exportMethod Summary
#' @keywords internal 
setMethod(
  "Summary", "mcmcestfix",
  function(x, ..., na.rm = FALSE) {
    dopt <- getOption("digits")
    obj <- x
    K <- obj@K
    rnames <- .rownames.Mcmcestfix(obj)
    cnames <- c("Estimates", "Std. Error")
    cat("\n")
    cat("Call: mcmcestimate\n")
    cat("\n")
    cat("Method: Gibbs Sampling with fixed indicators\n")
    cat("\n")
    cat(paste("Number of Iterations: ", obj@M, "\n", sep = ""))
    cat(paste("Number of Burnin Iterations: ", obj@burnin,
      "\n",
      sep = ""
    ))
    cat("\n")
    cat("Parameters:\n")
    cat("\n")
    cat(paste("Component Parameters: ",
      .parnames.Mcmcestfix(obj), "\n",
      sep = ""
    ))
    ## MAP ##
    cat("Maximum A Posterior (MAP)\n")
    parout <- .pars.map.Mcmcestfix(obj)
    rownames(parout) <- rnames
    colnames(parout) <- cnames
    print(parout)
    cat("\n")
    cat(paste("Log likelihood: ", sprintf("%.4f", obj@map$log), "\n", sep = ""))
    cat("---\n")
    ## BML ##
    cat("Bayesian Maximum Likelihood (BML)\n")
    parout <- .pars.bml.Mcmcestfix(obj)
    rownames(parout) <- rnames
    colnames(parout) <- cnames
    print(parout)
    cat("\n")
    cat(paste("Log likelihood: ", sprintf("%.4f", obj@bml$log), "\n", sep = ""))
    cat("---\n")
    ## IEAVG ##
    cat("Identified Ergodic Average (IEAVG)\n")
    parout <- .pars.ieavg.Mcmcestfix(obj)
    rownames(parout) <- rnames
    colnames(parout) <- cnames
    print(parout)
    cat("---\n")
    options(digits = dopt)
  }
)

## Getters ##
#' Getter method of `mcmcestfix` class.
#' 
#' Returns the `dist` slot.
#' 
#' @param object An `mcmcestfix` object.
#' @returns The `dist` slot of the `object`.
#' @exportMethod getDist
#' @keywords internal
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
#' f_est <- mcmcestimate(f_output)
#' # Get the slot.
#' getDist(f_est)
#' 
#' @seealso 
#' * [mcmcestind-class] for the corresponding class for models 
#'   with unknown indicators
#' * [mcmcestimate()] for calculating point estimates from MCMC samples
setMethod(
  "getDist", "mcmcestfix",
  function(object) {
    return(object@dist)
  }
)

#' Getter method of `mcmcestfix` class.
#' 
#' Returns the `K` slot.
#' 
#' @param object An `mcmcestfix` object.
#' @returns The `K` slot of the `object`.
#' @exportMethod getK
#' @keywords internal
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
#' f_est <- mcmcestimate(f_output)
#' # Get the slot.
#' getK(f_est)
#' 
#' @seealso 
#' * [mcmcestind-class] for the corresponding class for models 
#'   with unknown indicators
#' * [mcmcestimate()] for calculating point estimates from MCMC samples
setMethod(
  "getK", "mcmcestfix",
  function(object) {
    return(object@K)
  }
)

#' Getter method of `mcmcestfix` class.
#' 
#' Returns the `indicmod` slot.
#' 
#' @param object An `mcmcestfix` object.
#' @returns The `indicmod` slot of the `object`.
#' @exportMethod getIndicmod
#' @keywords internal
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
#' f_est <- mcmcestimate(f_output)
#' # Get the slot.
#' getIndicmod(f_est)
#' 
#' @seealso 
#' * [mcmcestind-class] for the corresponding class for models 
#'   with unknown indicators
#' * [mcmcestimate()] for calculating point estimates from MCMC samples
setMethod(
  "getIndicmod", "mcmcestfix",
  function(object) {
    return(object@indicmod)
  }
)

#' Getter method of `mcmcestfix` class.
#' 
#' Returns the `burnin` slot.
#' 
#' @param object An `mcmcestfix` object.
#' @returns The `burnin` slot of the `object`.
#' @exportMethod getBurnin
#' @keywords internal
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
#' f_est <- mcmcestimate(f_output)
#' # Get the slot.
#' getBurnin(f_est)
#' 
#' @seealso 
#' * [mcmcestind-class] for the corresponding class for models 
#'   with unknown indicators
#' * [mcmcestimate()] for calculating point estimates from MCMC samples
setMethod(
  "getBurnin", "mcmcestfix",
  function(object) {
    return(object@burnin)
  }
)

#' Getter method of `mcmcestfix` class.
#' 
#' Returns the `M` slot.
#' 
#' @param object An `mcmcestfix` object.
#' @returns The `M` slot of the `object`.
#' @exportMethod getM
#' @keywords internal
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
#' f_est <- mcmcestimate(f_output)
#' # Get the slot.
#' getM(f_est)
#' 
#' @seealso 
#' * [mcmcestind-class] for the corresponding class for models 
#'   with unknown indicators
#' * [mcmcestimate()] for calculating point estimates from MCMC samples
setMethod(
  "getM", "mcmcestfix",
  function(object) {
    return(object@M)
  }
)

#' Getter method of `mcmcestfix` class.
#' 
#' Returns the `ranperm` slot.
#' 
#' @param object An `mcmcestfix` object.
#' @returns The `ranperm` slot of the `object`.
#' @exportMethod getRanperm
#' @keywords internal
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
#' f_est <- mcmcestimate(f_output)
#' # Get the slot.
#' getRanperm(f_est)
#' 
#' @seealso 
#' * [mcmcestind-class] for the corresponding class for models 
#'   with unknown indicators
#' * [mcmcestimate()] for calculating point estimates from MCMC samples
setMethod(
  "getRanperm", "mcmcestfix",
  function(object) {
    return(object)
  }
)

#' Getter method of `mcmcestfix` class.
#' 
#' Returns the `relabel` slot.
#' 
#' @param object An `mcmcestfix` object.
#' @returns The `relabel` slot of the `object`.
#' @exportMethod getRelabel
#' @keywords internal
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
#' f_est <- mcmcestimate(f_output)
#' # Get the slot.
#' getRelabel(f_est)
#' 
#' @seealso 
#' * [mcmcestind-class] for the corresponding class for models 
#'   with unknown indicators
#' * [mcmcestimate()] for calculating point estimates from MCMC samples
setMethod(
  "getRelabel", "mcmcestfix",
  function(object) {
    return(object@relabel)
  }
)

#' Getter method of `mcmcestfix` class.
#' 
#' Returns the `map` slot.
#' 
#' @param object An `mcmcestfix` object.
#' @returns The `map` slot of the `object`.
#' @exportMethod getMap
#' @keywords internal
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
#' f_est <- mcmcestimate(f_output)
#' # Get the slot.
#' getMap(f_est)
#' 
#' @seealso 
#' * [mcmcestind-class] for the corresponding class for models 
#'   with unknown indicators
#' * [mcmcestimate()] for calculating point estimates from MCMC samples
setMethod(
  "getMap", "mcmcestfix",
  function(object) {
    return(object@map)
  }
)

#' Getter method of `mcmcestfix` class.
#' 
#' Returns the `bml` slot.
#' 
#' @param object An `mcmcestfix` object.
#' @returns The `bml` slot of the `object`.
#' @exportMethod getBml
#' @keywords internal
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
#' f_est <- mcmcestimate(f_output)
#' # Get the slot.
#' getBml(f_est)
#' 
#' @seealso 
#' * [mcmcestind-class] for the corresponding class for models 
#'   with unknown indicators
#' * [mcmcestimate()] for calculating point estimates from MCMC samples
setMethod(
  "getBml", "mcmcestfix",
  function(object) {
    return(object@bml)
  }
)

#' Getter method of `mcmcestfix` class.
#' 
#' Returns the `ieavg` slot.
#' 
#' @param object An `mcmcestfix` object.
#' @returns The `ieavg` slot of the `object`.
#' @exportMethod getIeavg
#' @keywords internal
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
#' f_est <- mcmcestimate(f_output)
#' # Get the slot.
#' getIeavg(f_est)
#' 
#' @seealso 
#' * [mcmcestind-class] for the corresponding class for models 
#'   with unknown indicators
#' * [mcmcestimate()] for calculating point estimates from MCMC samples
setMethod(
  "getIeavg", "mcmcestfix",
  function(object) {
    return(object@ieavg)
  }
)

#' Getter method of `mcmcestfix` class.
#' 
#' Returns the `sdpost` slot.
#' 
#' @param object An `mcmcestfix` object.
#' @returns The `sdpost` slot of the `object`.
#' @exportMethod getSdpost
#' @keywords internal
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
#' f_est <- mcmcestimate(f_output)
#' # Get the slot.
#' getIeavg(f_est)
#' 
#' @seealso 
#' * [mcmcestind-class] for the corresponding class for models 
#'   with unknown indicators
#' * [mcmcestimate()] for calculating point estimates from MCMC samples
setMethod(
  "getSdpost", "mcmcestfix",
  function(object) {
    return(object@sdpost)
  }
)

## No setters as users are not intended to manipulate
## this object

### Private functions.
### These functions are not exported.

### Summary

#' Summarize MAP estimates
#' 
#' @description 
#' For internal usage only. This function generates explicit summaries for the 
#' MAP estimates. 
#' 
#' Note that at this time advanced summaries are only available for Poisson 
#' mixture models.
#' 
#' @param obj An `mcmcestfix` object containing the parameter estimates.
#' @return A matrix with parameter estimates from the MAP. In addition the 
#'   standard deviations of the posterior density are presented.
#' @noRd 
#' 
#' @seealso 
#' * [Summary()] for the calling function
".pars.map.Mcmcestfix" <- function(obj) {
  if (obj@dist == "poisson") {
    .pars.map.poisson.Mcmcestfix(obj)
  }
}

#' Summarize MAP estimates form Poisson mixture models
#' 
#' @description 
#' For internal usage only. This function generates explicit summaries for the 
#' MAP estimates of Poisson mixture models.
#' 
#' @param obj An `mcmcestfix` object containing the parameter estimates.
#' @return A matrix with parameter estimates from the MAP. In addition the 
#'   standard deviations of the posterior density are presented.
#' @noRd 
#' 
#' @seealso 
#' * [Summary()] for the calling function
".pars.map.poisson.Mcmcestfix" <- function(obj) {
  parout <- matrix(0, nrow = obj@K, ncol = 2)
  for (k in seq(1, obj@K)) {
    parout[k, 1] <- obj@map$par$lambda[k]
    parout[k, 2] <- obj@sdpost$identified$par$lambda[k]
  }
  return(parout)
}

#' Summarize BML estimates
#' 
#' @description 
#' For internal usage only. This function generates explicit summaries for the 
#' BML estimates. 
#' 
#' Note that at this time advanced summaries are only available for Poisson 
#' mixture models.
#' 
#' @param obj An `mcmcestfix` object containing the parameter estimates.
#' @return A matrix with parameter estimates from the BML. In addition the 
#'   standard deviations of the posterior density are presented.
#' @noRd 
#' 
#' @seealso 
#' * [Summary()] for the calling function
".pars.bml.Mcmcestfix" <- function(obj) {
  if (obj@dist == "poisson") {
    .pars.bml.poisson.Mcmcestfix(obj)
  }
}

#' Summarize BML estimates for Poisson mixture models
#' 
#' @description 
#' For internal usage only. This function generates explicit summaries for the 
#' BML estimates. 
#' 
#' @param obj An `mcmcestfix` object containing the parameter estimates.
#' @return A matrix with parameter estimates from the BML. In addition the 
#'   standard deviations of the posterior density are presented.
#' @noRd 
#' 
#' @seealso 
#' * [Summary()] for the calling function
".pars.bml.poisson.Mcmcestfix" <- function(obj) {
  parout <- matrix(0, nrow = obj@K, ncol = 2)
  for (k in seq(1, obj@K)) {
    parout[k, 1] <- obj@bml$par$lambda[k]
    parout[k, 2] <- obj@sdpost$identified$par$lambda[k]
  }
  return(parout)
}

#' Summarize IEAVG estimates
#' 
#' @description 
#' For internal usage only. This function generates explicit summaries for the 
#' IEAVG estimates. 
#' 
#' Note that at this time advanced summaries are only available for Poisson 
#' mixture models.
#' 
#' @param obj An `mcmcestfix` object containing the parameter estimates.
#' @return A matrix with parameter estimates from the IEAVG. In addition the 
#'   standard deviations of the posterior density are presented.
#' @noRd 
#' 
#' @seealso 
#' * [Summary()] for the calling function
".pars.ieavg.Mcmcestfix" <- function(obj) {
  if (obj@dist == "poisson") {
    .pars.ieavg.poisson.Mcmcestfix(obj)
  }
}

#' Summarize IEAVG estimates for Poisson mixture models
#' 
#' @description 
#' For internal usage only. This function generates explicit summaries for the 
#' IEAVG estimates. 
#' 
#' @param obj An `mcmcestfix` object containing the parameter estimates.
#' @return A matrix with parameter estimates from the IEAVG. In addition the 
#'   standard deviations of the posterior density are presented.
#' @noRd 
#' 
#' @seealso 
#' * [Summary()] for the calling function
".pars.ieavg.poisson.Mcmcestfix" <- function(obj) {
  parout <- matrix(0, nrow = obj@K, ncol = 2)
  for (k in seq(1, obj@K)) {
    parout[k, 1] <- obj@ieavg$par$lambda[k]
    parout[k, 2] <- obj@sdpost$identified$par$lambda[k]
  }
  return(parout)
}

#' Create summary row names
#' 
#' @description 
#' For internal usage only. This function generates row names for the explicit 
#' summaries.
#' 
#' Note that at this time advanced summaries are only available for Poisson 
#' mixture models.
#' 
#' @param obj An `mcmcestfix` object containing the parameter estimates.
#' @return A vector with the row names for the advanced summary.
#' @noRd 
#' 
#' @seealso 
#' * [Summary()] for the calling function
".rownames.Mcmcestfix" <- function(obj) {
  if (obj@dist == "poisson") {
    .rownames.poisson.Mcmcestfix(obj)
  }
}

#' Create summary row names for Poisson mixture models
#' 
#' @description 
#' For internal usage only. This function generates row names for the explicit 
#' summaries.
#' 
#' @param obj An `mcmcestfix` object containing the parameter estimates.
#' @return A vector with the row names for the advanced summary over estimates 
#'   for a Poisson mixture model.
#' @noRd 
#' 
#' @seealso 
#' * [Summary()] for the calling function
".rownames.poisson.Mcmcestfix" <- function(obj) {
  rnames <- rep("", obj@K)
  for (k in seq(1, obj@K)) {
    rnames[k] <- paste("lambda ", k, sep = "")
  }
  return(rnames)
}

#' Create parameter names for components
#' 
#' @description 
#' For internal usage only. This function generates parameter names to be used 
#' in the advanced summary. 
#' 
#' Note that at this time advanced summaries are only available for Poisson 
#' mixture models.
#' 
#' @param obj An `mcmcestfix` object containing the parameter estimates.
#' @return A vector with the names of the component parameters to be used in 
#'   the advanced summary.
#' @noRd 
#' 
#' @seealso 
#' * [Summary()] for the calling function
".parnames.Mcmcestfix" <- function(obj) {
  if (obj@dist == "poisson") {
    parnames <- c("lambda")
  }
  return(parnames)
}