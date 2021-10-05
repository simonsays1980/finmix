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

#' Finmix `cdatamoments` class
#' 
#' Stores moments of an [fdata][fdata_class] object containing continuous data. 
#' The `fdata` object is stored in the parent [datamoments][datamoments_class] 
#' class.
#' 
#' @slot higher An array containing the four higher centralized moments of the 
#'   continuous data stored in the `fdata` object.
#' @slot skewness A vector storing the skewness of the continuous data in the 
#'   corresponding `fdata` object. 
#' @slot kurtosis A vector storing the kurtosis of the continuous data in the 
#'   corresponding `fdata` object. 
#' @slot corr A matrix containing the correlations between the data dimensions 
#'   in case of multivariate data (i.e. slot `r` in the `fdata` object is 
#'   larger than one).
#' @slot smoments A `csdatamoments` object, if the `fdata` object also holds 
#'   indicators. `NULL`, if no indicators are present in the `fdata` object. 
#' @exportClass cdatamoments
#' @name cdatamoments_class
#' @seealso 
#' * [datamoments][datamoments_class] for the parent class
#' * [ddatamoments][ddatamoments_class] for the corresponding class for 
#'   discrete data
#' * [csdatamoments][csdatamoments_class] for the contained class if indicators
#'   are present in the `fdata` object 
.cdatamoments <- setClass("cdatamoments",
  representation(
    higher = "array",
    skewness = "vector",
    kurtosis = "vector",
    corr = "matrix",
    smoments = "csdatamomentsOrNULL"
  ),
  contains = c("datamoments"),
  validity = function(object) {
    ## else: ok
    TRUE
  },
  prototype(
    higher = array(),
    skewness = vector(),
    kurtosis = vector(),
    corr = matrix(),
    smoments = .csdatamoments()
  )
)

#' Initializer of the `cdatamoments` class
#' 
#' @description
#' Only used implicitly. The initializer calls a function `generateMoments()` 
#' to generate in the initialization step the moments for a passed-in `fdata`
#' object.
#' 
#' @param .Object An object: see the "initialize Methods" section in 
#'   [initialize].
#' @param ... Arguments to specify properties of the new object, to be passed 
#'   to `initialize()`.
#' @param model A finmix `fdata` object containing the observations.
#' @noRd
#' 
#' @seealso 
#' * [Classes_Details] for details of class definitions, and 
#' * [setOldClass] for the relation to S3 classes
setMethod(
  "initialize", "cdatamoments",
  function(.Object, ..., value = fdata()) {
    .Object@fdata <- value
    if (hasS(value)) {
      .Object@smoments <- sdatamoments(value = value)
    } else {
      .Object@smoments <- NULL
    }
    generateMoments(.Object)
  }
)

#' Generate moments for continuous data.
#' 
#' @description 
#' Implicit method. Calling [generateMoments()] generates the moments of a
#' finite mixture with continuous data.
#' 
#' @param object An `cdatamoments` object. 
#' @return An `cdatamoments` object with calculated moments.
#' @noRd
setMethod(
  "generateMoments", "cdatamoments",
  function(object) {
    .generateCdatamoments(object)
  }
)

#' Shows a summary of a `cdatamoments` object.
#' 
#' Calling [show()] on a `cdatamoments` object gives an overview 
#' of the moments of a finit mixture with continuous data.
#' 
#' @param object A `cdatamoments` object.
#' @returns A console output listing the slots and summary information about
#'   each of them. 
#' @exportMethod show
#' @describeIn cdatamoments_class
setMethod(
  "show", "cdatamoments",
  function(object) {
    cat("Object 'datamoments'\n")
    cat(
      "     mean        : Vector of",
      length(object@mean), "\n"
    )
    cat(
      "     var         : Vector of",
      length(object@var), "\n"
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
    if (!all(is.na(object@corr))) {
      cat(
        "     corr        :",
        paste(dim(object@corr), collapse = "x"), "\n"
      )
    }
    if (hasS(object@fdata)) {
      cat(
        "     smoments    : Object of class",
        class(object@smoments), "\n"
      )
    }
    cat(
      "     fdata        : Object of class",
      class(object@fdata), "\n"
    )
  }
)

## Getters ##
#' Getter method of `cdatamoments` class.
#' 
#' Returns the `smoments` slot.
#' 
#' @param object An `cdatamoments` object.
#' @returns The `smoments` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' # Generate an exponential mixture model with two components.
#' f_model <- model("exponential", par = list(lambda = c(0.3, 0.7)), K = 2)
#' # Simulate data from the model.
#' f_data <- simulate(f_model)
#' # Calculate the mixture moments.
#' f_datamoms <- datamoments(f_data)
#' # Get the moments for the included indicators of the data. 
#' getSmoments(f_datamoms)
#' 
#' @seealso 
#' * [datamoments][datamoments_class] for the base class for model moments
#' * [datamoments()] for the constructor of the `datamoments` class family
setMethod(
  "getSmoments", "cdatamoments",
  function(object) {
    return(object@smoments)
  }
)

#' Getter method of `cdatamoments` class.
#' 
#' Returns the `higher` slot.
#' 
#' @param object An `cdatamoments` object.
#' @returns The `higher` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' # Generate an exponential mixture model with two components.
#' f_model <- model("exponential", par = list(lambda = c(0.3, 0.7)), K = 2)
#' # Simulate data from the model.
#' f_data <- simulate(f_model)
#' # Calculate the mixture moments.
#' f_datamoms <- datamoments(f_data)
#' # Use the getter.
#' getHigher(f_datamoms)
#' 
#' @seealso 
#' * [datamoments][datamoments_class] for the base class for model moments
#' * [datamoments()] for the constructor of the `datamoments` class family
setMethod(
  "getHigher", "cdatamoments",
  function(object) {
    return(object@higher)
  }
)

#' Getter method of `cdatamoments` class.
#' 
#' Returns the `skewness` slot.
#' 
#' @param object An `cdatamoments` object.
#' @returns The `skewness` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' # Generate an exponential mixture model with two components.
#' f_model <- model("exponential", par = list(lambda = c(0.3, 0.7)), K = 2)
#' # Simulate data from the model.
#' f_data <- simulate(f_model)
#' # Calculate the mixture moments.
#' f_datamoms <- datamoments(f_data)
#' # Use the getter.
#' getSkewness(f_datamoms)
#' 
#' @seealso 
#' * [datamoments][datamoments_class] for the base class for model moments
#' * [datamoments()] for the constructor of the `datamoments` class family
setMethod(
  "getSkewness", "cdatamoments",
  function(object) {
    return(object@skewness)
  }
)

#' Getter method of `cdatamoments` class.
#' 
#' Returns the `kurtosis` slot.
#' 
#' @param object An `cdatamoments` object.
#' @returns The `kurtosis` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' # Generate an exponential mixture model with two components.
#' f_model <- model("exponential", par = list(lambda = c(0.3, 0.7)), K = 2)
#' # Simulate data from the model.
#' f_data <- simulate(f_model)
#' # Calculate the mixture moments.
#' f_datamoms <- datamoments(f_data)
#' # Use the getter.
#' getKurtosis(f_datamoms)
#' 
#' @seealso 
#' * [datamoments][datamoments_class] for the base class for model moments
#' * [datamoments()] for the constructor of the `datamoments` class family
setMethod(
  "getKurtosis", "cdatamoments",
  function(object) {
    return(object@kurtosis)
  }
)

#' Getter method of `cdatamoments` class.
#' 
#' Returns the `corr` slot.
#' 
#' @param object An `cdatamoments` object.
#' @returns The `corr` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' # Generate an exponential mixture model with two components.
#' f_model <- model("exponential", par = list(lambda = c(0.3, 0.7)), K = 2)
#' # Simulate data from the model.
#' f_data <- simulate(f_model)
#' # Calculate the mixture moments.
#' f_datamoms <- datamoments(f_data)
#' # Use the getter.
#' getCorr(f_datamoms)
#' 
#' @seealso 
#' * [datamoments][datamoments_class] for the base class for model moments
#' * [datamoments()] for the constructor of the `datamoments` class family
setMethod(
  "getCorr", "cdatamoments",
  function(object) {
    return(object@corr)
  }
)
## Setters ##
## No setters as users should not manipulate a 'cdatamoments' object ##

### Private functions
### These function are not exported
#' Generate data moments for finite mixture data
#' 
#' @description 
#' Only called implicitly. generates all moments of finite mixture data in a 
#' `fdata` object.
#' 
#' @param object An `cdatamoments` object to contain all calculated
#'   moments. 
#' @returns An `cdatamoments` object containing all moments of the 
#'   inite mixture data.
#' @importFrom stats var cor
#' @noRd
".generateCdatamoments" <- function(object) {
  ## enforce column-wise ordering ##
  hasY(object@fdata, verbose = TRUE)
  datam <- getColY(object@fdata)
  ## Compute higher moments ##
  ## higher.moments is a r x L matrix (L = 4) ##
  means <- apply(datam, 2, mean, na.rm = TRUE)
  object@mean <- means
  object@var <- var(datam, na.rm = TRUE)
  d <- datam - rep(means, each = nrow(datam))
  momentsm <- array(0, dim = c(4, object@fdata@r))
  momentsm[2, ] <- apply(d^2, 2, mean, na.rm = TRUE)
  momentsm[3, ] <- apply(d^3, 2, mean, na.rm = TRUE)
  momentsm[4, ] <- apply(d^4, 2, mean, na.rm = TRUE)
  dimnames(momentsm) <- list(
    c("1st", "2nd", "3rd", "4th"),
    colnames(datam)
  )
  object@higher <- momentsm
  ## Compute skewness and kurtosis ##
  ## skewness and kurtosis are 1 x r vectors ##
  skewm <- momentsm[3, ] / momentsm[2, ]^1.5
  kurtm <- momentsm[4, ] / momentsm[2, ]^2
  names(skewm) <- colnames(datam)
  names(kurtm) <- colnames(datam)
  object@skewness <- skewm
  object@kurtosis <- kurtm
  ## Compute corr matrix in case of r > 1 ##
  ## corr is a r x r matrix ##
  if (object@fdata@r > 1) {
    object@corr <- cor(datam)
  } else {
    object@corr <- matrix()
  }
  return(object)
}
