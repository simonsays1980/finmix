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

#' Finmix `studentmodelmoments` class
#' 
#' @description
#' Defines a class that holds theoretical moments for a finite mixture of 
#' student distributions. Note that this class is not directly used, but 
#' indirectly when calling the `modelmoments` constructor [modelmoments()].
#' 
#' @slot B A numeric defining the between-group heterogeneity.
#' @slot W A numeric defining the within-group heterogeneity. 
#' @slot R A numeric defining the coefficient of determination.
#' @exportClass studentmodelmoments
#' @name studentmodelmoments
#' 
#' @seealso 
#' * \code{\link{modelmoments_class}} for the base class for model moments
#' * \code{\link{modelmoments}} for the constructor of `modelmoments` classes
.studentmodelmoments <- setClass("studentmodelmoments",
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

#' Initializer of the `studentmodelmoments` class
#' 
#' @description
#' Only used implicitly. The initializer calls a function `generateMoments()` 
#' to generate in the initialization step also the moments for a passed `model`
#' object.
#' 
#' @param .Object An object: see the "initialize Methods" section in 
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
  "initialize", "studentmodelmoments",
  function(.Object, ..., model) {
    .Object <- callNextMethod(.Object, ..., model = model)
    generateMoments(.Object)
  }
)

#' Generate moments for student mixture
#' 
#' @description 
#' Implicit method. Calling [generateMoments()] generates the moments of an
#' student mixture distribution.
#' 
#' @param object An `studentmodelmoments` object. 
#' @return An `studentmodelmoments` object with calculated moments.
#' @noRd
setMethod(
  "generateMoments", "studentmodelmoments",
  function(object) {
    .generateMomentsStudent(object)
  }
)

#' Shows a summary of an `studentmodelmoments` object.
#' 
#' Calling [show()] on an `studentmodelmoments` object gives an overview 
#' of the moments of an student finite mixture. 
#' 
#' @param object An `studentmodelmoments` object.
#' @returns A console output listing the slots and summary information about
#'   each of them. 
#' @noRd
#' @seealso
#' * [modelmoments-class] for the base class for model moments
#' * [modelmoments()] for the constructor of `modelmoments` classes
setMethod(
  "show", "studentmodelmoments",
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
      paste(dim(object@higher), collapse = "x"), "\n"
    )
    cat(
      "     skewness    : Vector of",
      length(object@skewness), "\n"
    )
    cat(
      "     kurtosis    : Vector of",
      length(object@kurtosis), "\n"
    )
    cat("     B           :", object@B, "\n")
    cat("     W           :", object@W, "\n")
    cat("     R           :", object@R, "\n")
    cat(
      "     model       : Object of class",
      class(object@model), "\n"
    )
  }
)

## Getters ##
#' Getter method of `studentmodelmoments` class.
#' 
#' Returns the `B` slot.
#' 
#' @param object An `studentmodelmoments` object.
#' @returns The `B` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' f_model         <- model("normal", weight = matrix(c(.3, .7), nrow = 1))
#' means           <- c(-2, 2)
#' sigmas          <- matrix(c(2, 4), nrow=1)
#' setPar(f_model) <- list(mu = means, sigma = sigmas, df = c(20, 40))
#' f_moments       <- modelmoments(f_model)
#' getB(f_moments)
#' 
#' @seealso 
#' * [modelmoments] for the base class for model moments
#' * [modelmoments()] for the constructor of the `modelmoments` class family
setMethod(
  "getB", "studentmodelmoments",
  function(object) {
    return(object@B)
  }
)

#' Getter method of `studentmodelmoments` class.
#' 
#' Returns the `W` slot.
#' 
#' @param object An `studentmodelmoments` object.
#' @returns The `W` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' f_model         <- model("normal", weight = matrix(c(.3, .7), nrow = 1))
#' means           <- c(-2, 2)
#' sigmas          <- matrix(c(2, 4), nrow=1)
#' setPar(f_model) <- list(mu = means, sigma = sigmas, df = c(20, 40))
#' f_moments       <- modelmoments(f_model)
#' getW(f_moments)
#' 
#' @seealso 
#' * [modelmoments] for the base class for model moments
#' * [modelmoments()] for the constructor of the `modelmoments` class family
setMethod(
  "getW", "studentmodelmoments",
  function(object) {
    return(object@W)
  }
)

#' Getter method of `studentmodelmoments` class.
#' 
#' Returns the `R` slot.
#' 
#' @param object An `studentmodelmoments` object.
#' @returns The `R` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' f_model         <- model("normal", weight = matrix(c(.3, .7), nrow = 1))
#' means           <- c(-2, 2)
#' sigmas          <- matrix(c(2, 4), nrow=1)
#' setPar(f_model) <- list(mu = means, sigma = sigmas, df = c(20, 40))
#' f_moments       <- modelmoments(f_model)
#' getR(f_moments)
#' 
#' @seealso 
#' * [modelmoments] for the base class for model moments
#' * [modelmoments()] for the constructor of the `modelmoments` class family
setMethod(
  "getR", "studentmodelmoments",
  function(object) {
    return(object@R)
  }
)

## No setters as users are not intended to manipulate ##
## this object ##

### Private functions
### These function are not exported
".generateMomentsStudent" <- function(object) {
  mu <- object@model@par$mu
  sigma <- object@model@par$sigma
  df <- object@model@par$df
  weight <- object@model@weight
  object@mean <- sum(weight * mu)
  object@higher <- .mixturemoments.student(
    object@model,
    4, object@mean
  )
  dimnames(object@higher) <- list(
    c("1st", "2nd", "3rd", "4th"),
    ""
  )
  object@var <- array(object@higher[2], dim = c(1, 1))
  object@skewness <- object@higher[3] / object@higher[2]^1.5
  object@kurtosis <- object@higher[4] / object@higher[2]^2
  object@B <- sum(weight * (mu - object@mean)^2)
  object@W <- sum(weight * sigma * df / (df - 2))
  object@R <- 1 - object@W / object@var[1]
  return(object)
}
