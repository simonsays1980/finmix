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
#' Finmix `exponentialmodelmoments` class
#' 
#' @description
#' Defines a class that holds modelmoments for a finite mixture of exponential
#' distributions. Note that this class is not directly used, but indirectly 
#' when calling the `modelmoments` constructor [modelmoments()].
#' 
#' @slot B A numeric defining the between-group heterogeneity.
#' @slot W A numeric defining the within-group heterogeneity. 
#' @slot R A numeric defining the coefficient of determination.
#' @exportClass exponentialmodelmoments
#' @name exponentialmodelmoments
#' 
#' @seealso 
#' * [modelmoments_class] for the base class for model moments
#' * [modelmoments()] for the constructor of `modelmoments` classes
.exponentialmodelmoments <- setClass("exponentialmodelmoments",
  representation(
    B = "numeric",
    W = "numeric",
    R = "numeric"
  ),
  contains = c("cmodelmoments"),
  validity = function(object) {
    ## else: OK
    TRUE
  },
  prototype(
    B = numeric(),
    W = numeric(),
    R = numeric()
  )
)

#' Initializer of the `exponentialmoments` class
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
  "initialize", "exponentialmodelmoments",
  function(.Object, ..., model) {
    .Object <- callNextMethod(.Object, ..., model = model)
    generateMoments(.Object)
  }
)

#' Generate moments for exponential mixture
#' 
#' @description 
#' Implicit method. Calling [generateMoments()] generates the moments of an
#' exponential mixture distribution.
#' 
#' @param object An `exponentialmodelmoments` object. 
#' @return An `exponentialmodelmoments` object with calculated moments.
#' @noRd
setMethod(
  "generateMoments", "exponentialmodelmoments",
  function(object) {
    .generateMomentsExponential(object)
  }
)

#' Shows a summary of an `exponentialmodelmoments` object.
#' 
#' Calling [show()] on an `exponentialmodelmoments` object gives an overview 
#' of the moments of an exponential finite mixture. 
#' 
#' @param object An `exponentialmodelmoments` object.
#' @returns A console output listing the slots and summary information about
#'   each of them. 
#' @describeIn exponentialmodelmoments
setMethod(
  "show", "exponentialmodelmoments",
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
      "     higher      :",
      paste(dim(object@higher), collapse = "x"),
      "\n"
    )
    cat("     skewness    :", object@skewness, "\n")
    cat("     kurtosis    :", object@kurtosis, "\n")
    cat("     B           :", object@B, "\n")
    cat("     W           :", object@W, "\n")
    cat("     R           :", object@R, "\n")
    cat(
      "     model       : Object of class",
      class(object@model), "\n"
    )
  }
)

#' Getter method of `exponentialmodelmoments` class.
#' 
#' Returns the `B` slot.
#' 
#' @param object An `exponentialmodelmoments` object.
#' @returns The `B` slot of the `object`.
#' @describeIn modelmoments Getter method for slot `B`
#' 
#' @examples 
#' f_model <- model("exponential", par=list(lambda=c(0.3, 0.1)), 
#'                  weight=matrix(c(0.3, 0.7), nrow=1))
#' f_moments <- modelmoments(f_model)
#' getB(f_moments)
#' 
#' @seealso 
#' * [modelmoments] for the base class for model moments
#' * [modelmoments()] for the constructor of the `modelmoments` class family
setMethod(
  "getB", "exponentialmodelmoments",
  function(object) {
    return(object@B)
  }
)

#' Getter method of `exponentialmodelmoments` class.
#' 
#' Returns the `W` slot.
#' 
#' @param object An `exponentialmodelmoments` object.
#' @returns The `W` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' f_model <- model("exponential", par=list(lambda=c(0.3, 0.1)), 
#'                  weight=matrix(c(0.3, 0.7), nrow=1))
#' f_moments <- modelmoments(f_model)
#' getW(f_moments)
#' 
#' @seealso 
#' * [modelmoments] for the base class for model moments
#' * [modelmoments()] for the constructor of the `modelmoments` class family
setMethod(
  "getW", "exponentialmodelmoments",
  function(object) {
    return(object@W)
  }
)

#' Getter method of `exponentialmodelmoments` class.
#' 
#' Returns the `R` slot.
#' 
#' @param object An `exponentialmodelmoments` object.
#' @returns The `R` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' f_model <- model("exponential", par=list(lambda=c(0.3, 0.1)), 
#'                  weight=matrix(c(0.3, 0.7), nrow=1))
#' f_moments <- modelmoments(f_model)
#' getR(f_moments)
#' 
#' @seealso 
#' * [modelmoments] for the base class for model moments
#' * [modelmoments()] for the constructor of the `modelmoments` class family
setMethod(
  "getR", "exponentialmodelmoments",
  function(object) {
    return(object@R)
  }
)

## No setters as users are not intended to manipulate ##
## this object ##

### Private functions
### These functions are not exported
#' Generate model moments for an exponential mixture
#' 
#' @description 
#' Only called implicitly. generates all moments of an exponential mixture 
#' distribution.
#' 
#' @param object An `exponentialmodelmoments` object to contain all calculated
#'   moments. 
#' @returns An `exponentialmodelmoments` object containing all moments of the 
#'   exponential mixture distributions.
#' @noRd
".generateMomentsExponential" <- function(object) {
  lambda <- object@model@par$lambda
  weight <- object@model@weight
  object@mean <- sum(weight * 1 / lambda)
  highm <- .mixturemoments.exponential(object@model, 4, object@mean)
  dimnames(highm) <- list(c("1st", "2nd", "3rd", "4th"), "")
  object@higher <- highm
  object@var <- array(object@higher[2], dim = c(1, 1))
  object@skewness <- object@higher[3] / object@higher[2]^1.5
  object@kurtosis <- object@higher[4] / object@higher[2]^2
  object@W <- sum(weight * 1 / lambda^2)
  object@B <- sum(weight * (1 / lambda - object@mean)^2)
  object@R <- 1 - object@W / object@var[1]
  return(object)
}
