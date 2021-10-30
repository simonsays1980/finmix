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
#' Finmix `normultmodelmoments` class
#' 
#' @description
#' Defines a class that holds modelmoments for a finite mixture of normult
#' distributions. Note that this class is not directly used, but indirectly 
#' when calling the `modelmoments` constructor [modelmoments()].
#' 
#' @slot B A numeric defining the between-group heterogeneity.
#' @slot W A numeric defining the within-group heterogeneity. 
#' @slot R A numeric defining the coefficient of determination.
#' @exportClass normultmodelmoments
#' @name normultmodelmoments-class
#' 
#' @seealso 
#' * [modelmoments-class] for the base class for model moments
#' * [modelmoments()] for the constructor of `modelmoments` classes
.normultmodelmoments <- setClass("normultmodelmoments",
  representation(
    B = "array",
    W = "array",
    Rdet = "numeric",
    Rtr = "numeric",
    corr = "array"
  ),
  contains = c("cmodelmoments"),
  validity = function(object) {
    ## else: OK
    TRUE
  },
  prototype(
    B    = array(),
    W    = array(),
    Rdet = numeric(),
    Rtr  = numeric(),
    corr = array()
  )
)

#' Initializer of the `normultmoments` class
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
  "initialize", "normultmodelmoments",
  function(.Object, ..., model) {
    .Object <- callNextMethod(.Object, ..., model = model)
    generateMoments(.Object)
  }
)

#' Generate moments for normult mixture
#' 
#' @description 
#' Implicit method. Calling [generateMoments()] generates the moments of an
#' normult mixture distribution.
#' 
#' @param object An `normultmodelmoments` object. 
#' @return An `normultmodelmoments` object with calculated moments.
#' @noRd
setMethod(
  "generateMoments", "normultmodelmoments",
  function(object) {
    .generateMomentsNormult(object)
  }
)

#' Shows a summary of an `normultmodelmoments` object.
#' 
#' Calling [show()] on an `normultmodelmoments` object gives an overview 
#' of the moments of an normult finite mixture. 
#' 
#' @param object An `normultmodelmoments` object.
#' @returns A console output listing the slots and summary information about
#'   each of them. 
#' @exportMethod show
#' @keywords internal
#' @seealso 
#' * [modelmoments-class] for the base class for model moments
#' * [modelmoments()] for the constructor of `modelmoments` classes
setMethod(
  "show", "normultmodelmoments",
  function(object) {
    cat("Object 'modelmoments'\n")
    cat(
      "     mean        : Vector of",
      length(object@mean), "\n"
    )
    cat(
      "     var         :",
      paste(dim(object@var), collapse = "x"),
      "\n"
    )
    cat(
      "     higher      :",
      paste(dim(object@higher), collapse = "x"),
      "\n"
    )
    cat(
      "     skewness    : Vector of",
      length(object@skewness), "\n"
    )
    cat(
      "     kurtosis    : Vector of",
      length(object@kurtosis), "\n"
    )
    cat(
      "     B           :",
      paste(dim(object@B), collapse = "x"), "\n"
    )
    cat(
      "     W           :",
      paste(dim(object@W), collapse = "x"), "\n"
    )
    cat("     Rdet        :", object@Rdet, "\n")
    cat("     Rtr         :", object@Rtr, "\n")
    cat(
      "     corr        :",
      paste(dim(object@corr), collapse = "x"),
      "\n"
    )
    cat(
      "     model       : Object of class",
      class(object@model), "\n"
    )
  }
)

## Getters ##
#' Getter method of `normultmodelmoments` class.
#' 
#' @description 
#' Returns the `B` slot.
#' 
#' @param object An `normultmodelmoments` object.
#' @return The `B` slot of the `object`.
#' @exportMethod getB
#' @keywords internal
#' 
#' @examples 
#' f_model         <- model("normult", weight = matrix(c(.3, .7), nrow = 1))
#' means           <- matrix(c(-2, -2, 2, 2),nrow = 2)
#' covar           <- matrix(c(1, 1.2, 1.2, 4), nrow = 2)
#' sigmas          <- array(c(covar, 2*covar), dim = c(2, 2, 2))
#' setPar(f_model) <- list(mu = means, sigma = sigmas)
#' f_moments       <- modelmoments(f_model)
#' getB(f_moments)
#' 
#' @seealso 
#' * [modelmoments] for the base class for model moments
#' * [modelmoments()] for the constructor of the `modelmoments` class family
setMethod(
  "getB", "normultmodelmoments",
  function(object) {
    return(object@B)
  }
)

#' Getter method of `normultmodelmoments` class.
#' 
#' Returns the `W` slot.
#' 
#' @param object An `normultmodelmoments` object.
#' @returns The `W` slot of the `object`.
#' @exportMethod getW
#' @keywords internal
#' 
#' @examples 
#' f_model         <- model("normult", weight = matrix(c(.3, .7), nrow = 1))
#' means           <- matrix(c(-2, -2, 2, 2),nrow = 2)
#' covar           <- matrix(c(1, 1.2, 1.2, 4), nrow = 2)
#' sigmas          <- array(c(covar, 2*covar), dim = c(2, 2, 2))
#' setPar(f_model) <- list(mu = means, sigma = sigmas)
#' f_moments       <- modelmoments(f_model)
#' getW(f_moments)
#' 
#' @seealso 
#' * [modelmoments] for the base class for model moments
#' * [modelmoments()] for the constructor of the `modelmoments` class family
setMethod(
  "getW", "normultmodelmoments",
  function(object) {
    return(object@W)
  }
)

#' Getter method of `normultmodelmoments` class.
#' 
#' Returns the `Rdet` slot.
#' 
#' @param object An `normultmodelmoments` object.
#' @returns The `Rdet` slot of the `object`.
#' @exportMethod getRdet
#' @keywords internal
#' 
#' @examples 
#' f_model         <- model("normult", weight = matrix(c(.3, .7), nrow = 1))
#' means           <- matrix(c(-2, -2, 2, 2),nrow = 2)
#' covar           <- matrix(c(1, 1.2, 1.2, 4), nrow = 2)
#' sigmas          <- array(c(covar, 2*covar), dim = c(2, 2, 2))
#' setPar(f_model) <- list(mu = means, sigma = sigmas)
#' f_moments       <- modelmoments(f_model)
#' getRdet(f_moments)
#' 
#' @seealso 
#' * [modelmoments] for the base class for model moments
#' * [modelmoments()] for the constructor of the `modelmoments` class family
setMethod(
  "getRdet", "normultmodelmoments",
  function(object) {
    return(object@Rdet)
  }
)

#' Getter method of `normultmodelmoments` class.
#' 
#' Returns the `Rtr` slot.
#' 
#' @param object An `normultmodelmoments` object.
#' @returns The `Rtr` slot of the `object`.
#' @exportMethod getRtr
#' @keywords internal
#' 
#' @examples 
#' f_model         <- model("normult", weight = matrix(c(.3, .7), nrow = 1))
#' means           <- matrix(c(-2, -2, 2, 2),nrow = 2)
#' covar           <- matrix(c(1, 1.2, 1.2, 4), nrow = 2)
#' sigmas          <- array(c(covar, 2*covar), dim = c(2, 2, 2))
#' setPar(f_model) <- list(mu = means, sigma = sigmas)
#' getRtr(f_moments)
#' 
#' @seealso 
#' * [modelmoments] for the base class for model moments
#' * [modelmoments()] for the constructor of the `modelmoments` class family
setMethod(
  "getRtr", "normultmodelmoments",
  function(object) {
    return(object@B)
  }
)

#' Getter method of `normultmodelmoments` class.
#' 
#' Returns the `Corr` slot.
#' 
#' @param object An `normultmodelmoments` object.
#' @returns The `Corr` slot of the `object`.
#' @exportMethod getCorr
#' @keywords internal
#' 
#' @examples 
#' f_model         <- model("normult", weight = matrix(c(.3, .7), nrow = 1))
#' means           <- matrix(c(-2, -2, 2, 2),nrow = 2)
#' covar           <- matrix(c(1, 1.2, 1.2, 4), nrow = 2)
#' sigmas          <- array(c(covar, 2*covar), dim = c(2, 2, 2))
#' setPar(f_model) <- list(mu = means, sigma = sigmas)
#' f_moments       <- modelmoments(f_model) 
#' getCorr(f_moments)
#' 
#' @seealso 
#' * [modelmoments] for the base class for model moments
#' * [modelmoments()] for the constructor of the `modelmoments` class family
setMethod(
  "getCorr", "normultmodelmoments",
  function(object) {
    return(object@corr)
  }
)

## No setters as users are not intended to manipulate ##
## this object ##

### private functions
### these function are not exported
#' Generate model moments for an normult mixture
#' 
#' @description 
#' Only called implicitly. generates all moments of an normult mixture 
#' distribution.
#' 
#' @param object An `normultmodelmoments` object to contain all calculated
#'   moments. 
#' @returns An `normultmodelmoments` object containing all moments of the 
#'   normult mixture distributions.
#' @noRd
".generateMomentsNormult" <- function(object) {
  mu <- object@model@par$mu
  sigma <- object@model@par$sigma
  weight <- object@model@weight
  names <- rep("", object@model@r)
  for (i in seq(1, object@model@r)) {
    names[i] <- paste("r=", i, sep = "")
  }
  object@mean <- apply(apply(mu, 1, "*", weight),
    2, sum,
    na.rm = TRUE
  )
  object@W <- apply(sweep(sigma, MARGIN = 1, weight, "*"),
    c(1, 2), sum,
    na.rm = TRUE
  )
  object@var <- object@W + apply(
    apply(mu, 2, tcrossprod, mu),
    1, "*", weight
  )
  object@var <- object@var - object@mean %*% t(object@mean)
  diffm <- mu - object@mean
  object@B <- apply(
    apply(diffm, 1, tcrossprod, diffm),
    1, "*", weight
  )
  cd <- diag(1 / diag(object@var)^.5)
  object@corr <- cd %*% object@var %*% cd
  object@Rtr <- 1 - sum(diag(object@W)) / sum(diag(object@var))
  object@Rdet <- 1 - det(object@W) / det(object@var)
  highm <- array(0, dim = c(4, object@model@r))
  for (i in seq(1, object@model@r)) {
    marmodel <- mixturemar(object@model, i)
    highm[, i] <- t(.mixturemoments.normal(
      marmodel,
      4,
      object@mean[i]
    ))
  }
  names(object@mean) <- names
  colnames(object@var) <- names
  rownames(object@var) <- names
  colnames(object@B) <- names
  rownames(object@B) <- names
  colnames(object@W) <- names
  rownames(object@W) <- names
  colnames(object@corr) <- names
  rownames(object@corr) <- names
  object@higher <- highm
  dimnames(object@higher) <- list(c("1st", "2nd", "3rd", "4th"), names)
  object@skewness <- object@higher[3, ] / object@higher[2, ]^1.5
  object@kurtosis <- object@higher[4, ] / object@higher[2, ]^2
  return(object)
}
