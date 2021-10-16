# Copyright (C) 2013 Lars Simon Zehnder
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

#' Finmix `groupmoments` class
#' 
#' Stores moments for finite mixture component distributions. These are only
#' available, if the data contains in addition to observations also indicators 
#' defining to which component a certain observation belongs. These indicators 
#' are stored in an [fdata-class] object in the slot `S`. 
#' 
#' @slot NK An array containing the group sizes for each component.
#' @slot mean A matrix containing the group averages for each component.
#' @slot WK An array containing the within-group variability. For multivariate 
#'   data this is an array of dimension `K x r x r` and for univariate 
#'   data this is simply an array of dimension `1 x K`.
#' @slot var An array containing the within-group (co)variance. For multivariate 
#'   data this is an array of dimension `K x r x r` and for univariate 
#'   data this is simply an array of dimension `1 x K`.
#' @slot fdata An [fdata-class] object containing the data. 
#' @exportClass groupmoments
#' @rdname groupmoments-class
#' @keywords internal
#' @seealso 
#' * [groupmoments()] for the class constructor
#' * [datamoments-class] for the base class for data moments
#' * [datamoments()] for the constructor of any object of the `datamoments` 
#'   class family
.groupmoments <- setClass("groupmoments",
  representation(
    NK = "array",
    mean = "matrix",
    WK = "array",
    var = "array",
    fdata = "fdata"
  ),
  validity = function(object) {
    ## else: ok
    TRUE
  },
  prototype(
    NK = array(),
    mean = matrix(),
    WK = array(),
    var = array(),
    fdata = fdata()
  )
)

#' Finmix `groupmoments` class constructor
#' 
#' @description
#' Calling [groupmoments()] creates an object holding various 
#' component-specific moments. These moments can only constructed if the 
#' [fdata-class] object contains in addition to observations also 
#' indicators defining from which component a certain observation stems.
#' 
#' @param value An `fdata` object containing observations in slot `y` and 
#'   indicators in slot `S`.
#' @return A `groupmoments` object containing component-specific moments of the 
#'   `fdata` object.
#' @export
#' @name groupmoments
#' 
#' @examples 
#' # Define a mixture model with exponential components.
#' f_model <- model("exponential", par = list(lambda = c(0.3, 0.7)), K = 2)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Create group moments of the data.
#' groupmoments(f_data)
#' 
#' @seealso 
#' * [fdata-class] for the `fdata` class definition
#' * [groupmoments-class] for the definition of the `groupmoments` 
#'   class
#' * [datamoments-class] for the base class for data moments
#' * [datamoments()] for the constructor of any object of the `datamoments` 
#'   class family
"groupmoments" <- function(value = fdata()) {
  hasY(value, verbose = TRUE)
  hasS(value, verbose = TRUE)
  .groupmoments(value = value)
}

## initializes by immediately calling method ##
## 'generateMoments' ##
#' Initializer of the `groupmoments` class
#' 
#' @description
#' Only used implicitly. The initializer calls a function `generateMoments()` 
#' object. to generate in the initialization step the moments for a passed-in 
#' `fdata` object.
#' 
#' @param .Object An object: see the "initialize Methods" section in 
#'   [initialize].
#' @param ... Arguments to specify properties of the new object, to be passed 
#'   to `initialize()`.
#' @param model A finmix [fdata-class] object containing the observations.
#' @noRd
#' 
#' @seealso 
#' * [Classes_Details] for details of class definitions, and 
#' * [setOldClass] for the relation to S3 classes
setMethod(
  "initialize", "groupmoments",
  function(.Object, ..., value) {
    .Object@fdata <- value
    generateMoments(.Object)
  }
)

#' Generate moments
#' 
#' @description 
#' Implicit method. Calling [generateMoments()] generates the moments of a
#' finite mixture with continuous data.
#' 
#' @param object A `groupmoments` object. 
#' @return An `groupmoments` object with calculated moments.
#' @noRd
setMethod(
  "generateMoments", "groupmoments",
  function(object) {
    .generateGroupMoments(object)
  }
)

## R usual 'show' function ##
#' Shows a summary of a `groupmoments` object.
#' 
#' Calling [show()] on a `groupmoments` object gives an overview 
#' of the moments of a finit mixture with continuous data.
#' 
#' @param object A `groupmoments` object.
#' @returns A console output listing the slots and summary information about
#'   each of them. 
#' @exportMethod show
#' @noRd
setMethod(
  "show", "groupmoments",
  function(object) {
    cat("Object 'groupmoments'\n")
    cat(
      "     NK          : Vector of",
      length(object@NK), "\n"
    )
    cat(
      "     mean        :",
      paste(dim(object@mean), collapse = "x"), "\n"
    )
    cat(
      "     WK          :",
      paste(dim(object@WK), collapse = "x"), "\n"
    )
    cat(
      "     var         :",
      paste(dim(object@var), collapse = "x"), "\n"
    )
    cat(
      "     fdata       : Object of class",
      class(object@fdata), "\n"
    )
  }
)

## R usual Getters ##
#' Getter method of `groupmoments` class.
#' 
#' Returns the `NK` slot.
#' 
#' @param object An `groupmoments` object.
#' @returns The `NK` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' # Generate a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 0.7)), K = 2)
#' # Simulate data from the model.
#' f_data <- simulate(f_model)
#' # Calculate the mixture moments.
#' f_gmoments <- groupmoments(f_data)
#' # Get the moments for the included indicators of the data. 
#' getNK(f_gmoments)
#' 
#' @seealso 
#' * [groupmoments][groupmments_class] for the definition of the `groupmoments` 
#'   class
#' * [groupmoments()] for the class constructor
setMethod(
  "getNK", "groupmoments",
  function(object) {
    return(object@NK)
  }
)

#' Getter method of `groupmoments` class.
#' 
#' Returns the `mean` slot.
#' 
#' @param object An `groupmoments` object.
#' @returns The `mean` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' # Generate a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 0.7)), K = 2)
#' # Simulate data from the model.
#' f_data <- simulate(f_model)
#' # Calculate the mixture moments.
#' f_gmoments <- groupmoments(f_data)
#' # Get the moments for the included indicators of the data. 
#' getMean(f_gmoments)
#' 
#' @seealso 
#' * [groupmoments][groupmments_class] for the definition of the `groupmoments` 
#'   class
#' * [groupmoments()] for the class constructor
setMethod(
  "getMean", "groupmoments",
  function(object) {
    return(object@mean)
  }
)

#' Getter method of `groupmoments` class.
#' 
#' Returns the `WK` slot.
#' 
#' @param object An `groupmoments` object.
#' @returns The `WK` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' # Generate a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 0.7)), K = 2)
#' # Simulate data from the model.
#' f_data <- simulate(f_model)
#' # Calculate the mixture moments.
#' f_gmoments <- groupmoments(f_data)
#' # Get the moments for the included indicators of the data. 
#' getWK(f_gmoments)
#' 
#' @seealso 
#' * [groupmoments][groupmments_class] for the definition of the `groupmoments` 
#'   class
#' * [groupmoments()] for the class constructor
setMethod(
  "getWK", "groupmoments",
  function(object) {
    return(object@WK)
  }
)

#' Getter method of `groupmoments` class.
#' 
#' Returns the `Var` slot.
#' 
#' @param object An `groupmoments` object.
#' @returns The `Var` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' # Generate a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 0.7)), K = 2)
#' # Simulate data from the model.
#' f_data <- simulate(f_model)
#' # Calculate the mixture moments.
#' f_gmoments <- groupmoments(f_data)
#' # Get the moments for the included indicators of the data. 
#' getVar(f_gmoments)
#' 
#' @seealso 
#' * [groupmoments][groupmments_class] for the definition of the `groupmoments` 
#'   class
#' * [groupmoments()] for the class constructor
setMethod(
  "getVar", "groupmoments",
  function(object) {
    return(object@var)
  }
)

setMethod(
  "getFdata", "groupmoments",
  function(object) {
    return(object@fdata)
  }
)
## No setters as user are not intended to manipulate this  ##
## object ##

### Private functions
### These functions are not exported
#' Generate data moments for finite mixture data
#' 
#' @description 
#' Only called implicitly. generates all moments of finite mixture data in a 
#' `fdata` object.
#' 
#' @param object A `groupmoments` object to contain all calculated
#'   moments. 
#' @returns A `groupmoments` object containing all moments of the 
#'   finite mixture data.
#' @noRd
".generateGroupMoments" <- function(object) {
  if (!hasS(object@fdata)) {
    return(object)
  }

  ## Compute group sizes ##
  ## enforce column-wise ordering ##

  datam <- getColY(object@fdata)
  classm <- getColS(object@fdata)
  ## Calculate group sizes and group means ##
  ## 'NK' is an 1 x K vector ##
  ## 'groupmean' is an r x K matrix ##
  level.set <- as.numeric(levels(factor(classm)))
  K <- length(level.set)
  r <- ncol(datam)
  comp <- matrix(rep(classm, K), ncol = K) == matrix(seq(1, K),
    nrow = nrow(datam),
    ncol = K,
    byrow = TRUE
  )
  names <- rep("", K)
  for (k in seq(1, K)) {
    names[k] <- paste("k=", k, sep = "")
  }
  object@NK <- as.array(apply(comp, 2, sum))
  dimnames(object@NK) <- list(names)
  gmeans <- matrix(NA, nrow = r, ncol = K)
  for (i in seq(1, r)) {
    gmeans[i, ] <- (t(datam[, i]) %*% comp) / t(object@NK)
  }
  colnames(gmeans) <- names
  rownames(gmeans) <- colnames(datam)
  object@mean <- gmeans
  wkm <- array(NA, dim = c(r, r, K))
  varm <- array(NA, dim = c(r, r, K))
  for (k in seq(1, K)) {
    group.demeaned <- (datam - rep(gmeans[, k], each = nrow(datam))) * comp[, k]
    wkm[, , k] <- t(group.demeaned) %*% group.demeaned
    varm[, , k] <- wkm[, , k] / object@NK[k]
  }
  dimnames(wkm) <- list(colnames(datam), colnames(datam), names)
  dimnames(varm) <- list(colnames(datam), colnames(datam), names)
  object@WK <- wkm
  object@var <- varm
  return(object)
}
