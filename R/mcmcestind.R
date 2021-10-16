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
#' Note that this class inherits almost all of its slots from the `mcmcestfix` 
#' class, the corresponding class for fixed indicators.
#' 
#' @slot eavg A named list containing the estimates of the ergodic average. The 
#'   element `par` is a list and contains the component parameter estimates and 
#'   `weight` contains the weight estimates. The difference between the EAVG 
#'   and the IEAVG is that the IEAVG is based on re-labeled samples.
#' @exportClass mcmcestind
#' @rdname mcmcestind-class
#' @keywords internal
#' 
#' @seealso
#' * [mcmcestfix-class] for the parent class with fixed indicators
#' * [mcmcestimate()] to calculate point estimates
.mcmcestind <- setClass("mcmcestind",
  representation(eavg = "list"),
  contains = c("mcmcestfix"),
  validity = function(object) {
    ## else: OK
    TRUE
  },
  prototype(eavg = list())
)

#' Finmix `mcmcest` class union
#' 
#' @description 
#' This class union includes all classes that define objects for storing the 
#' parameter estimates and is used to dispatch methods for `mcmcest` objects.
#' 
#' @exportClass mcmcest
#' @name mcmcest-class
#' @noRd
setClassUnion(
  "mcmcest",
  c(
    "mcmcestfix",
    "mcmcestind"
  )
)

#' Shows a summary of an `mcmcestind` object.
#' 
#' Calling [show()] on an `mcmcestind` object gives an overview 
#' of the `mcmcestind` object.
#' 
#' @param object An `mcmcestind` object.
#' @returns A console output listing the slots and summary information about
#'   each of them. 
#' @exportMethod show
#' @noRd
setMethod(
  "show", "mcmcestind",
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
      "     eavg        : List of",
      length(object@eavg), "\n"
    )
    cat(
      "     sdpost      : List of",
      length(object@sdpost), "\n"
    )
  }
)

# TODO: The Std. Error is the same for both components.
#' Shows an advanced summary of an `mcmcestind` object.
#' 
#' Calling [Summary()] on an `mcmcestind` object gives an advanced overview 
#' of the `mcmcestind` object.
#' 
#' Note, this method is so far only implemented for mixtures of Poisson 
#' distributions.
#' 
#' @param object An `mcmcestind` object.
#' @returns A console output listing the formatted slots and summary 
#'   information about each of them. 
#' @exportMethod Summary
setMethod(
  "Summary", "mcmcestind",
  function(x, ..., na.rm = FALSE) {
    dopt <- getOption("digits")
    options(digits = 4)
    obj <- x
    K <- obj@K
    rnames <- .rownames.Mcmcestind(obj)
    cnames <- c("Estimates", "Std. Error")
    cat("\n")
    cat("Call: mcmcestimate\n")
    cat("\n")
    if (obj@ranperm) {
      cat("Method: Random Permutation Gibbs Sampling\n")
    } else {
      cat("Method: Gibbs Sampling\n")
    }
    cat("\n")
    cat(paste("Number of Iterations: ", obj@M, "\n", sep = ""))
    cat(paste("Number of Burnin Iterations: ", obj@burnin,
      "\n",
      sep = ""
    ))
    cat(paste("Relabeling algorithm used: ", obj@relabel, "\n",
      sep = ""
    ))
    cat("\n")
    cat("Parameters:\n")
    cat("\n")
    cat(paste("Component Parameters: ",
      .parnames.Mcmcestfix(obj), "\n",
      sep = ""
    ))
    cat("Weights: eta\n")
    ## MAP ##
    cat("Maximum A Posterior (MAP)\n")
    parout <- .pars.map.Mcmcestind(obj)
    rownames(parout) <- rnames
    colnames(parout) <- cnames
    print(parout)
    cat("\n")
    cat(paste("Log likelihood: ", sprintf("%.4f", obj@map$log), "\n", sep = ""))
    cat("---\n")
    ## BML ##
    cat("Bayesian Maximum Likelihood (BML)\n")
    parout <- .pars.bml.Mcmcestind(obj)
    rownames(parout) <- rnames
    colnames(parout) <- cnames
    print(parout)
    cat("\n")
    cat(paste("Log likelihood: ", sprintf("%.4f", obj@bml$log), "\n", sep = ""))
    cat("---\n")
    ## IEAVG ##
    cat("Identified Ergodic Average (IEAVG)\n")
    parout <- .pars.ieavg.Mcmcestind(obj)
    rownames(parout) <- rnames
    colnames(parout) <- cnames
    print(parout)
    cat("---\n")
    ## EAVG ##
    cat("Ergodic Average (EAVG)\n")
    parout <- .pars.eavg.Mcmcestind(obj)
    rownames(parout) <- rnames
    colnames(parout) <- cnames
    print(parout)
    cat("---\n")
    options(digits = dopt)
  }
)

## Getters ##
#' Getter method of `mcmcestind` class.
#' 
#' Returns the `eavg` slot.
#' 
#' @param object An `mcmcestind` object.
#' @returns The `eavg` slot of the `object`.
#' @exportMethod getEavg
#' @noRd
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
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' f_est <- mcmcestimate(f_output)
#' # Get the slot.
#' getEavg(f_output)
#' 
#' @seealso 
#' * [mcmcestfix-class] for the parent class with fixed indicators
#' * [mcmcestimate()] for calculating point estimates from MCMC samples
setMethod(
  "getEavg", "mcmcestind",
  function(object) {
    return(object@eavg)
  }
)

## No setters as users are not intended to manipulate
## this object.

### Private functions.
### These functions are not exported.

### Summary

#' Summarize MAP estimates
#' 
#' @description 
#' For internal usage only. This function generates explicit summaries for the 
#' MAP estimates of models with unknown indicators. 
#' 
#' Note that at this time advanced summaries are only available for Poisson 
#' mixture models.
#' 
#' @param obj An `mcmcestind` object containing the parameter estimates.
#' @return A matrix with parameter estimates from the MAP.
#' @noRd 
#' 
#' @seealso 
#' * [Summary()] for the calling function
".pars.map.Mcmcestind" <- function(obj) {
  if (obj@dist == "poisson") {
    .pars.map.poisson.Mcmcestind(obj)
  }
}

#' Summarize MAP estimates form Poisson mixture models
#' 
#' @description 
#' For internal usage only. This function generates explicit summaries for the 
#' MAP estimates of Poisson mixture models when indicators are unknown.
#' 
#' @param obj An `mcmcestind` object containing the parameter estimates.
#' @return A matrix with parameter estimates from the MAP. In addition the 
#'   standard deviations of the posterior density are presented.
#' @noRd 
#' 
#' @seealso 
#' * [Summary()] for the calling function
".pars.map.poisson.Mcmcestind" <- function(obj) {
  K <- obj@K
  parout <- matrix(0, nrow = 2 * K, ncol = 2)
  for (k in seq(1, K)) {
    parout[k, 1] <- obj@map$par$lambda[k]
    parout[k, 2] <- obj@sdpost$identified$par$lambda[k]
  }
  for (k in seq(1, K)) {
    parout[k + K, 1] <- obj@map$weight[k]
    parout[k + K, 2] <- obj@sdpost$identified$weight[k]
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
#' @param obj An `mcmcestind` object containing the parameter estimates.
#' @return A matrix with parameter estimates from the BML. In addition the 
#'   standard deviations of the posterior density are presented.
#' @noRd 
#' 
#' @seealso 
#' * [Summary()] for the calling function
".pars.bml.Mcmcestind" <- function(obj) {
  if (obj@dist == "poisson") {
    .pars.bml.poisson.Mcmcestind(obj)
  }
}

#' Summarize BML estimates for Poisson mixture models
#' 
#' @description 
#' For internal usage only. This function generates explicit summaries for the 
#' BML estimates. 
#' 
#' @param obj An `mcmcestind` object containing the parameter estimates.
#' @return A matrix with parameter estimates from the BML. In addition the 
#'   standard deviations of the posterior density are presented.
#' @noRd 
#' 
#' @seealso 
#' * [Summary()] for the calling function
".pars.bml.poisson.Mcmcestind" <- function(obj) {
  K <- obj@K
  parout <- matrix(0, nrow = 2 * K, ncol = 2)
  for (k in seq(1, K)) {
    parout[k, 1] <- obj@bml$par$lambda[k]
    parout[k, 2] <- obj@sdpost$identified$par$lambda[k]
  }
  for (k in seq(1, K)) {
    parout[k + K, 1] <- obj@bml$weight[k]
    parout[k + K, 2] <- obj@sdpost$identified$weight[k]
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
#' @param obj An `mcmcestind` object containing the parameter estimates.
#' @return A matrix with parameter estimates from the IEAVG. In addition the 
#'   standard deviations of the posterior density are presented.
#' @noRd 
#' 
#' @seealso 
#' * [Summary()] for the calling function
".pars.ieavg.Mcmcestind" <- function(obj) {
  if (obj@dist == "poisson") {
    .pars.ieavg.poisson.Mcmcestind(obj)
  }
}

#' Summarize IEAVG estimates for Poisson mixture models
#' 
#' @description 
#' For internal usage only. This function generates explicit summaries for the 
#' IEAVG estimates. 
#' 
#' @param obj An `mcmcestind` object containing the parameter estimates.
#' @return A matrix with parameter estimates from the IEAVG. In addition the 
#'   standard deviations of the posterior density are presented.
#' @noRd 
#' 
#' @seealso 
#' * [Summary()] for the calling function
".pars.ieavg.poisson.Mcmcestind" <- function(obj) {
  K <- obj@K
  parout <- matrix(0, nrow = 2 * K, ncol = 2)
  for (k in seq(1, K)) {
    parout[k, 1] <- obj@ieavg$par$lambda[k]
    parout[k, 2] <- obj@sdpost$identified$par$lambda[k]
  }
  for (k in seq(1, K)) {
    parout[k + K, 1] <- obj@ieavg$weight[k]
    parout[k + K, 2] <- obj@sdpost$identified$weight[k]
  }
  return(parout)
}

#' Summarize IEAVG estimates
#' 
#' @description 
#' For internal usage only. This function generates explicit summaries for the 
#' EAVG estimates. The difference between the EAVG and the IEAVG is that the 
#' IEAVG is based on re-labeled samples.
#' 
#' Note that at this time advanced summaries are only available for Poisson 
#' mixture models.
#' 
#' @param obj An `mcmcestind` object containing the parameter estimates.
#' @return A matrix with parameter estimates from the EAVG. In addition the 
#'   standard deviations of the posterior density are presented.
#' @noRd 
#' 
#' @seealso 
#' * [Summary()] for the calling function
".pars.eavg.Mcmcestind" <- function(obj) {
  if (obj@dist == "poisson") {
    .pars.eavg.poisson.Mcmcestind(obj)
  }
}

#' Summarize EAVG estimates for Poisson mixture models
#' 
#' @description 
#' For internal usage only. This function generates explicit summaries for the 
#' EAVG estimates. 
#' 
#' @param obj An `mcmcestind` object containing the parameter estimates.
#' @return A matrix with parameter estimates from the EAVG. In addition the 
#'   standard deviations of the posterior density are presented.
#' @noRd 
#' 
#' @seealso 
#' * [Summary()] for the calling function
".pars.eavg.poisson.Mcmcestind" <- function(obj) {
  K <- obj@K
  parout <- matrix(0, nrow = 2 * K, ncol = 2)
  for (k in seq(1, K)) {
    parout[k, 1] <- obj@eavg$par$lambda[k]
    parout[k, 2] <- obj@sdpost$unidentified$par$lambda[k]
  }
  for (k in seq(1, K)) {
    parout[k + K, 1] <- obj@eavg$weight[k]
    parout[k + K, 2] <- obj@sdpost$unidentified$weight[k]
  }

  return(parout)
}

#' Create row names for summary
#' 
#' @description 
#' For internal usage only. This function generates row names for the explicit 
#' summaries.
#' 
#' Note that at this time advanced summaries are only available for Poisson 
#' mixture models.
#' 
#' @param obj An `mcmcestind` object containing the parameter estimates.
#' @return A vector with the row names for the advanced summary.
#' @noRd 
#' 
#' @seealso 
#' * [Summary()] for the calling function
".rownames.Mcmcestind" <- function(obj) {
  if (obj@dist == "poisson") {
    .rownames.poisson.Mcmcestind(obj)
  }
}

#' Create summary row names for Poisson mixture models
#' 
#' @description 
#' For internal usage only. This function generates row names for the explicit 
#' summaries.
#' 
#' @param obj An `mcmcestind` object containing the parameter estimates.
#' @return A vector with the row names for the advanced summary over estimates 
#'   for a Poisson mixture model.
#' @noRd 
#' 
#' @seealso 
#' * [Summary()] for the calling function
".rownames.poisson.Mcmcestind" <- function(obj) {
  rnames <- rep("", 2 * obj@K)
  for (k in seq(1, obj@K)) {
    rnames[k] <- paste("lambda ", k, sep = "")
  }
  for (k in seq(obj@K + 1, 2 * obj@K)) {
    rnames[k] <- paste("eta ", k - obj@K, sep = "")
  }
  return(rnames)
}