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

#' Finmix `binomialmodelmoments` class
#' 
#' @description
#' Defines a class that holds modelmoments for a finite mixture of Binomial
#' distributions. Note that this class is not directly used, but indirectly 
#' when calling the `modelmoments` constructor [modelmoments()].
#' 
#' This is a class that directly inherits from the `dmodelmoments` class. 
#' @import methods
#' @exportClass binomialmodelmoments
#' @name binomialmodelmoments-class
#' @keywords internal 
#' @seealso 
#' * [modelmoments-class] for the base class for model moments
#' * [modelmoments()] for the constructor of `modelmoments` classes
#' * [dmodelmoments-class] class for the parent class
.binomialmodelmoments <- setClass("binomialmodelmoments",
  representation(extrabinvar = "numeric"),
  contains = c("dmodelmoments"),
  validity = function(object) {
    ## else: OK
    TRUE
  },
  prototype(extrabinvar = numeric())
)

#' Initializer of the `binomialmoments` class
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
  "initialize", "binomialmodelmoments",
  function(.Object, ..., model) {
    .Object <- callNextMethod(.Object, ..., model = model)
    generateMoments(.Object)
  }
)

#' Generate moments for binomial mixture
#' 
#' @description 
#' Implicit method. Calling [generateMoments()] generates the moments of an
#' binomial mixture distribution.
#' 
#' @param object An `binomialmodelmoments` object. 
#' @return An `binomialmodelmoments` object with calculated moments.
#' @noRd
setMethod(
  "generateMoments", "binomialmodelmoments",
  function(object) {
    .generateMomentsBinomial(object)
  }
)

#' Shows a summary of an `binomialmodelmoments` object.
#' 
#' Calling [show()] on an `binomialmodelmoments` object gives an overview 
#' of the moments of an binomial finite mixture. 
#' 
#' @param object An `binomialmodelmoments` object.
#' @returns A console output listing the slots and summary information about
#'   each of them. 
#' @exportMethod show
#' @noRd
#' @seealso 
#' * [modelmoments()] for the mutual constructor for all modelmoments
#' * [binomialmodelmoments-class] for the class definition 
setMethod(
  "show", "binomialmodelmoments",
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
      "     extrabinvar :", object@extrabinvar,
      "\n"
    )
    cat(
      "     model       : Object of class",
      class(object@model), "\n"
    )
  }
)

## Getters ##
#' Getter method of `binomialmodelmoments` class.
#' 
#' Returns the `extrabinvar` slot.
#' 
#' @param object An `binomialmodelmoments` object.
#' @returns The `extrabinvar` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' f_model <- model("binomial", par=list(p=c(0.3, 0.5)), 
#'                  weight=matrix(c(0.3, 0.7), nrow=1))
#' f_moments <- modelmoments(f_model)
#' getExtrabinvar(f_moments)
#' 
#' @seealso 
#' * \code{\link{modelmoments-class}} for the base class for model moments
#' * \code{\link{modelmoments}} for the constructor of the `modelmoments` class family
setMethod(
  "getExtrabinvar", "binomialmodelmoments",
  function(object) {
    return(object@extrabinvar)
  }
)

## No setters as users are not intended to manipulate ##
## this object ##

### Private functions
### These function are not exported
#' Generates theoretical moments for a binomial mixture
#' 
#' @description 
#' Calling [.genwerateMomentsBinomial()] generates theoretical model moments 
#' for the binomial model defined in the `model` object. Next to the general 
#' mixture moments available to any mixture model, the binomial moments also
#' include the extra-binomial variation `extrabinvar` 
#' (see Fruehwirth-Schnatter (2006)) and the number of expected zeros `zero`.
#' 
#' @param object A `binomialmodelmoments` object with correpsonding `model` 
#'   object. Note that, if the `model` object has repetitions in slot `T` with 
#'   dimension larger than one only the first repetition is used for theoretical
#'   moments. 
#' @returns A `modelmoments` object containing the theoreitcal moments of the 
#'   binomial mixture defined in the `model` object. 
#' @noRd
#' 
#' @seealso 
#' * [dmodelmoments-class] for the class definition of `dmodelmoments` 
#' * [modelmoments()] for the constructor calling this function 
".generateMomentsBinomial" <- function(object) {
  p <- object@model@par$p
  T <- object@model@T[1]
  weight <- object@model@weight
  object@mean <- sum(weight * p)
  object@var <- array(sum(weight * (T * p - object@mean)^2)
  + sum(weight * T * p * (1 - p)), dim = c(1, 1))
  factm <- array(NA, dim = c(4, 1))
  factm[1] <- object@mean
  for (i in seq(2, 4)) {
    if (T >= i) {
      factm[i] <- sum(weight * factorial(T) / factorial(T - i) * p^i)
    } else {
      factm[i] <- NaN
    }
  }
  dimnames(factm) <- list(c("1st", "2nd", "3rd", "4th"), "")
  object@factorial <- factm
  if (object@model@K > 1) {
    object@over <- object@var[1] - object@mean
  } else {
    object@over <- 0
  }
  object@zero <- sum(weight * (1 - p)^T)
  object@extrabinvar <- object@mean * (1 - object@mean / T)
  return(object)
}
