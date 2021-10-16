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
#' Finmix `poissonmodelmoments` class
#' 
#' @description
#' Defines a class that holds modelmoments for a finite mixture of poisson
#' distributions. Note that this class is not directly used, but indirectly 
#' when calling the `modelmoments` constructor [modelmoments()].
#' 
#' @exportClass poissonmodelmoments
#' @rdname poissonmodelmoments-class
#' @keywords internal
#' 
#' @seealso 
#' * [dmodelmoments-class] for the parent class 
#' * [modelmoments-class] for the base class for model moments
#' * [modelmoments()] for the constructor of `modelmoments` classes
.poissonmodelmoments <- setClass("poissonmodelmoments",
  contains = c("dmodelmoments"),
  validity = function(object) {
    ## else: OK
    TRUE
  }
)

#' Initializer of the `poissonmoments` class
#' 
#' @description
#' Only used implicitly. The initializer calls a function `generateMoments()` 
#' to generate in the initialization step also the moments for a passed `model`
#' object.
#' 
#' @param .Object An object_ see the "initialize Methods" section in 
#'   [initialize].
#' @param ... Arguments to specify properties of the new object, to be passed 
#'   to `initialize()`.
#' @param model A finmix `model` object containing the definition of the 
#'   finite mixture distribution.
#' @noRd
#' 
#' @seealso 
#' * [Classes_Details] for details of class definitions, and 
#' * [setOldClass] for the relation to S3 classes
setMethod(
  "initialize", "poissonmodelmoments",
  function(.Object, ..., model) {
    .Object <- callNextMethod(.Object, ..., model = model)
    generateMoments(.Object)
  }
)

#' Generate moments for poisson mixture
#' 
#' @description 
#' Implicit method. Calling [generateMoments()] generates the moments of an
#' poisson mixture distribution.
#' 
#' @param object An `poissonmodelmoments` object. 
#' @return An `poissonmodelmoments` object with calculated moments.
#' @noRd
setMethod(
  "generateMoments", "poissonmodelmoments",
  function(object) {
    .generateMomentsPoisson(object)
  }
)

#' Shows a summary of an `poissonmodelmoments` object.
#' 
#' Calling [show()] on an `poissonmodelmoments` object gives an overview 
#' of the moments of an poisson finite mixture. 
#' 
#' @param object An `poissonmodelmoments` object.
#' @returns A console output listing the slots and summary information about
#'   each of them. 
#' @exportMethod show
#' @noRd
#' @seealso
#' * [modelmoments-class] for the base class for model moments
#' * [modelmoments()] for the constructor of `modelmoments` classes
setMethod(
  "show", "poissonmodelmoments",
  function(object) {
    cat("Object 'modelmoments'\n")
    cat(
      "     mean        : Vector of",
      length(object@mean), "\n"
    )
    cat(
      "     var         :",
      paste(dim(object@var), collapse = "x"), "\n"
    )
    cat(
      "     factorial   :",
      paste(dim(object@factorial), collapse = "x"),
      "\n"
    )
    cat("     over        :", object@over, "\n")
    cat("     zero        :", object@zero, "\n")
    cat(
      "     model       : Object of class",
      class(object@model), "\n"
    )
  }
)

## No Setters as users are not intended to manipulate
## this object ##

### Private functions
### These functions are not exported
#' Generate model moments for an poisson mixture
#' 
#' @description 
#' Only called implicitly. generates all moments of an poisson mixture 
#' distribution.
#' 
#' @param object An `poissonmodelmoments` object to contain all calculated
#'   moments. 
#' @returns An `poissonmodelmoments` object containing all moments of the 
#'   poisson mixture distributions.
#' @noRd
".generateMomentsPoisson" <- function(object) {
  hasPar(object@model, verbose = TRUE)
  K <- object@model@K
  lambda <- object@model@par$lambda
  fact.names <- list(c("1st", "2nd", "3rd", "4th"), "")
  if (K == 1) {
    object@mean <- lambda
    object@var <- as.matrix(lambda)
    object@over <- 0
    factm <- array(NA, dim = c(4, 1))
    for (i in seq(1, 4)) {
      factm[i] <- lambda^i
    }
    dimnames(factm) <- fact.names
    object@factorial <- factm
    object@zero <- exp((-1) * lambda)
  } else {
    hasWeight(object@model, verbose = TRUE)
    weight <- object@model@weight
    object@mean <- sum(weight * lambda)
    object@var <- array(sum(weight * lambda * (lambda + 1))
    - object@mean^2, dim = c(1, 1))
    object@over <- object@var[1] - object@mean
    factm <- array(NA, dim = c(4, 1))
    for (i in seq(1, 4)) {
      factm[i] <- sum(weight * lambda^i)
    }
    dimnames(factm) <- fact.names
    object@factorial <- factm
    object@zero <- sum(weight * exp((-1) * lambda))
  }
  return(object)
}
