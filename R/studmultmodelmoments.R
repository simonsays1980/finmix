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
#' Finmix `studmultmodelmoments` class
#' 
#' @description
#' Defines a class that holds modelmoments for a finite mixture of studmult
#' distributions. Note that this class is not directly used, but indirectly 
#' when calling the `modelmoments` constructor [modelmoments()].
#' 
#' @slot B A numeric defining the between-group heterogeneity.
#' @slot W A numeric defining the within-group heterogeneity. 
#' @slot Rdet A numeric defining the coefficient of determination based on the
#'   determinant of the covariance matrix. 
#' @slot Rtr A numeric defining the coefficient of determination based on the 
#'   trace of the covariance matrix.
#' @slot corr A `matrix` storing the correlation matrix.
#' @exportClass studmultmodelmoments
#' @name studmultmodelmoments-class
#' 
#' @seealso 
#' * [modelmoments-class] for the base class for model moments
#' * [modelmoments()] for the constructor of `modelmoments` classes
.studmultmodelmoments <- setClass("studmultmodelmoments",
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
    B = array(),
    W = array(),
    Rdet = numeric(),
    Rtr = numeric(),
    corr = array()
  )
)

#' Initializer of the `studmultmoments` class
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
#' @keywords internal
#' 
#' @seealso 
#' * [Classes_Details] for details of class definitions, and 
#' * [setOldClass] for the relation to S3 classes
setMethod(
  "initialize", "studmultmodelmoments",
  function(.Object, ..., model) {
    .Object <- callNextMethod(.Object, ..., model = model)
    generateMoments(.Object)
  }
)

#' Generate moments for studmult mixture
#' 
#' @description 
#' Implicit method. Calling [generateMoments()] generates the moments of an
#' studmult mixture distribution.
#' 
#' @param object An `studmultmodelmoments` object. 
#' @return An `studmultmodelmoments` object with calculated moments.
#' @keywords internal
setMethod(
  "generateMoments", "studmultmodelmoments",
  function(object) {
    .generateMomentsStudmult(object)
  }
)

#' Shows a summary of an `studmultmodelmoments` object.
#' 
#' Calling [show()] on an `studmultmodelmoments` object gives an overview 
#' of the moments of an studmult finite mixture. 
#' 
#' @param object An `studmultmodelmoments` object.
#' @returns A console output listing the slots and summary information about
#'   each of them. 
#' @exportMethod show
#' @keywords internal
setMethod(
  "show", "studmultmodelmoments",
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
      paste(dim(object@corr), collapse = "x"), "\n"
    )
    cat(
      "     model       : Object of class",
      class(object@model), "\n"
    )
  }
)

## Getters ##
#' Getter method of `studmultmodelmoments` class.
#' 
#' @description 
#' Returns the `B` slot.
#' 
#' @param object An `studmultmodelmoments` object.
#' @return The `B` slot of the `object`.
#' @keywords internal
#' @exportMethod getB
#' 
#' @examples 
#' f_model         <- model("studmult", weight = matrix(c(.3, .7), nrow = 1))
#' means           <- matrix(c(-2, -2, 2, 2),nrow = 2)
#' covar           <- matrix(c(1, 1.2, 1.2, 4), nrow = 2)
#' sigmas          <- array(c(covar, 2*covar), dim = c(2, 2, 2))
#' setPar(f_model) <- list(mu = means, sigma = sigmas, df = c(10, 20))
#' f_moments       <- modelmoments(f_model)
#' getB(f_moments)
#' 
#' @seealso 
#' * [modelmoments] for the base class for model moments
#' * [modelmoments()] for the constructor of the `modelmoments` class family
setMethod(
  "getB", "studmultmodelmoments",
  function(object) {
    return(object@B)
  }
)

#' Getter method of `studmultmodelmoments` class.
#' 
#' Returns the `W` slot.
#' 
#' @param object An `studmultmodelmoments` object.
#' @returns The `W` slot of the `object`.
#' @exportMethod getW
#' @keywords internal
#' 
#' @examples 
#' f_model         <- model("studmult", weight = matrix(c(.3, .7), nrow = 1))
#' means           <- matrix(c(-2, -2, 2, 2),nrow = 2)
#' covar           <- matrix(c(1, 1.2, 1.2, 4), nrow = 2)
#' sigmas          <- array(c(covar, 2*covar), dim = c(2, 2, 2))
#' setPar(f_model) <- list(mu = means, sigma = sigmas, df = c(10, 20))
#' f_moments       <- modelmoments(f_model)
#' getW(f_moments)
#' 
#' @seealso 
#' * [modelmoments] for the base class for model moments
#' * [modelmoments()] for the constructor of the `modelmoments` class family
setMethod(
  "getW", "studmultmodelmoments",
  function(object) {
    return(object@W)
  }
)

#' Getter method of `studmultmodelmoments` class.
#' 
#' Returns the `Rdet` slot.
#' 
#' @param object An `studmultmodelmoments` object.
#' @returns The `Rdet` slot of the `object`.
#' @exportMethod getRdet
#' @keywords internal
#' 
#' @examples 
#' f_model         <- model("studmult", weight = matrix(c(.3, .7), nrow = 1))
#' means           <- matrix(c(-2, -2, 2, 2),nrow = 2)
#' covar           <- matrix(c(1, 1.2, 1.2, 4), nrow = 2)
#' sigmas          <- array(c(covar, 2*covar), dim = c(2, 2, 2))
#' setPar(f_model) <- list(mu = means, sigma = sigmas, df = c(10, 20))
#' f_moments       <- modelmoments(f_model)
#' getRdet(f_moments)
#' 
#' @seealso 
#' * [modelmoments] for the base class for model moments
#' * [modelmoments()] for the constructor of the `modelmoments` class family
setMethod(
  "getRdet", "studmultmodelmoments",
  function(object) {
    return(object@Rdet)
  }
)

#' Getter method of `studmultmodelmoments` class.
#' 
#' Returns the `Rtr` slot.
#' 
#' @param object An `studmultmodelmoments` object.
#' @returns The `Rtr` slot of the `object`.
#' @exportMethod getRtr
#' @keywords internal
#' 
#' @examples 
#' f_model         <- model("studmult", weight = matrix(c(.3, .7), nrow = 1))
#' means           <- matrix(c(-2, -2, 2, 2),nrow = 2)
#' covar           <- matrix(c(1, 1.2, 1.2, 4), nrow = 2)
#' sigmas          <- array(c(covar, 2*covar), dim = c(2, 2, 2))
#' setPar(f_model) <- list(mu = means, sigma = sigmas, df = c(10, 20))
#' f_moments       <- modelmoments(f_model)
#' getRtr(f_moments)
#' 
#' @seealso 
#' * [modelmoments] for the base class for model moments
#' * [modelmoments()] for the constructor of the `modelmoments` class family
setMethod(
  "getRtr", "studmultmodelmoments",
  function(object) {
    return(object@Rtr)
  }
)

#' Getter method of `studmultmodelmoments` class.
#' 
#' Returns the `Corr` slot.
#' 
#' @param object An `studmultmodelmoments` object.
#' @returns The `Corr` slot of the `object`.
#' @exportMethod getCorr
#' @keywords internal
#' 
#' @examples 
#' f_model         <- model("studmult", weight = matrix(c(.3, .7), nrow = 1))
#' means           <- matrix(c(-2, -2, 2, 2),nrow = 2)
#' covar           <- matrix(c(1, 1.2, 1.2, 4), nrow = 2)
#' sigmas          <- array(c(covar, 2*covar), dim = c(2, 2, 2))
#' setPar(f_model) <- list(mu = means, sigma = sigmas, df = c(10, 20))
#' f_moments       <- modelmoments(f_model)
#' getCorr(f_moments)
#' 
#' @seealso 
#' * [modelmoments] for the base class for model moments
#' * [modelmoments()] for the constructor of the `modelmoments` class family
setMethod(
  "getCorr", "studmultmodelmoments",
  function(object) {
    return(object@corr)
  }
)

## No setters as users are not intended to manipulate ##
## this object ##

### Private functions
### These function are not exported
#' Generate model moments for an studmult mixture
#' 
#' @description 
#' Only called implicitly. generates all moments of an studmult mixture 
#' distribution.
#' 
#' @param object An `studmultmodelmoments` object to contain all calculated
#'   moments. 
#' @returns An `studmultmodelmoments` object containing all moments of the 
#'   studmult mixture distributions.
#' @noRd
".generateMomentsStudmult" <- function(object) {
  mu <- object@model@par$mu
  sigma <- object@model@par$sigma
  df <- object@model@par$df
  weight <- object@model@weight
  names <- rep("", object@model@r)
  for (i in seq(1, object@model@r)) {
    names[i] <- paste("r=", i, sep = "")
  }
  object@mean <- apply(apply(mu, 1, "*", weight),
    2, sum,
    na.rm = TRUE
  )
  if (all(df > 2)) {
    object@W <- apply(sweep(sigma,
      MARGIN = 3,
      weight * df / (df - 2), "*"
    ),
    c(1, 2), sum,
    na.rm = TRUE
    )
    object@var <- object@W + apply(
      apply(
        mu, 2,
        tcrossprod, mu
      ),
      1, "*",
      weight
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
  } else {
    r <- object@model@r
    object@W <- array(NaN, dim = c(r, r))
    object@var <- array(NaN, dim = c(r, r))
    object@B <- array(NaN, dim = c(r, r))
    object@Rdet <- NaN
    object@Rtr <- NaN
    object@corr <- array(NaN, dim = c(r, r))
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
  highm <- array(0, dim = c(4, object@model@r))
  dimnames(highm) <- list(c("1st", "2nd", "3rd", "4th"), names)
  for (i in seq(1, object@model@r)) {
    marmodel <- mixturemar(object@model, i)
    highm[, i] <- .mixturemoments.student(
      marmodel, 4,
      object@mean[i]
    )
  }
  object@higher <- highm
  object@skewness <- object@higher[3, ] / object@higher[2, ]^1.5
  object@kurtosis <- object@higher[4, ] / object@higher[2, ]^2
  return(object)
}
