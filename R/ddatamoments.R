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
#' Finmix `ddatamoments` class
#' 
#' Stores moments of an [fdata-class] object containing discrete data. 
#' The `fdata` object is stored in the parent [datamoments-class] 
#' class.
#' 
#' @slot factorial An array containing the first four factorial moments of the 
#'   discrete data stored in the `fdata` object.
#' @slot over A vector storing the overdispersion of the discrete data in the 
#'   corresponding `fdata` object. 
#' @slot zero A vector storing the fractions of zeros in the observed data. <
#' @slot smoments An `sdatamoments` object, if the `fdata` object also holds 
#'   indicators. `NULL`, if no indicators are present in the `fdata` object. 
#' @exportClass ddatamoments
#' @rdname ddatamoments-class
#' @keywords internal
#' @seealso 
#' * [datamoments-class] for the parent class
#' * [ddatamoments-class] for the corresponding class for 
#'   continuous data
#' * [sdatamoments-class] for the contained class if indicators
#'   are present in the `fdata` object 
.ddatamoments <- setClass("ddatamoments",
  representation(
    factorial = "array",
    over = "vector",
    zero = "vector",
    smoments = "sdatamomentsOrNULL"
  ),
  contains = c("datamoments"),
  validity = function(object) {
    ## else: ok
    TRUE
  },
  prototype(
    factorial = array(),
    over = vector(),
    zero = vector(),
    smoments = .sdatamoments()
  )
)

#' Initializer of the `ddatamoments` class
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
  "initialize", "ddatamoments",
  function(.Object, ..., value = fdata()) {
    .Object@fdata <- value
    if (hasS(value)) {
      .Object@smoments <- sdatamoments(value)
    } else {
      .Object@smoments <- NULL
    }
    generateMoments(.Object)
  }
)

## Generic set in 'groupmoments.R' ##
#' Generate moments for continuous data.
#' 
#' @description 
#' Implicit method. Calling [generateMoments()] generates the moments of a
#' finite mixture with continuous data.
#' 
#' @param object An `ddatamoments` object. 
#' @return An `ddatamoments` object with calculated moments.
#' @noRd
setMethod(
  "generateMoments", "ddatamoments",
  function(object) {
    .generateDdatamoments(object)
  }
)

#' Shows a summary of a `ddatamoments` object.
#' 
#' Calling [show()] on a `ddatamoments` object gives an overview 
#' of the moments of a finit mixture with continuous data.
#' 
#' @param object A `ddatamoments` object.
#' @returns A console output listing the slots and summary information about
#'   each of them. 
#' @exportMethod show
#' @noRd
#' @seealso 
#' * [datamoments-class] for the parent class definition
#' * [datamoments()] for the mutual constructor of all datamoments classes
setMethod(
  "show", "ddatamoments",
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
      "     factorial   :",
      paste(dim(object@factorial), collapse = "x"), "\n"
    )
    cat(
      "     over        : Vector of",
      length(object@over), "\n"
    )
    cat(
      "     zero        : Vector of",
      length(object@zero), "\n"
    )
    if (hasS(object@fdata)) {
      cat(
        "     smoments    : Object of class",
        class(object@smoments), "\n"
      )
    }
    cat(
      "     fdata       : Object of class",
      class(object@fdata), "\n"
    )
  }
)

## Getters ##
#' Getter method of `ddatamoments` class.
#' 
#' Returns the `smoments` slot.
#' 
#' @param object An `ddatamoments` object.
#' @returns The `smoments` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' # Generate a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 0.7)), K = 2)
#' # Simulate data from the model.
#' f_data <- simulate(f_model)
#' # Calculate the mixture moments.
#' f_datamoms <- datamoments(f_data)
#' # Get the moments for the included indicators of the data. 
#' getSmoments(f_datamoms)
#' 
#' @seealso 
#' * [datamoments-class] for the base class for model moments
#' * [datamoments()] for the constructor of the `datamoments` class family
setMethod(
  "getSmoments", "ddatamoments",
  function(object) {
    return(object@smoments)
  }
)

#' Getter method of `ddatamoments` class.
#' 
#' Returns the `smoments` slot.
#' 
#' @param object An `ddatamoments` object.
#' @returns The `smoments` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' # Generate a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 0.7)), K = 2)
#' # Simulate data from the model.
#' f_data <- simulate(f_model)
#' # Calculate the mixture moments.
#' f_datamoms <- datamoments(f_data)
#' # Get the moments for the included indicators of the data. 
#' getFactorial(f_datamoms)
#' 
#' @seealso 
#' * [datamoments-class] for the base class for model moments
#' * [datamoments()] for the constructor of the `datamoments` class family
setMethod(
  "getFactorial", "ddatamoments",
  function(object) {
    return(object@factorial)
  }
)

#' Getter method of `ddatamoments` class.
#' 
#' Returns the `smoments` slot.
#' 
#' @param object An `ddatamoments` object.
#' @returns The `smoments` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' # Generate a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 0.7)), K = 2)
#' # Simulate data from the model.
#' f_data <- simulate(f_model)
#' # Calculate the mixture moments.
#' f_datamoms <- datamoments(f_data)
#' # Get the moments for the included indicators of the data. 
#' getOver(f_datamoms)
#' 
#' @seealso 
#' * [datamoments-class] for the base class for model moments
#' * [datamoments()] for the constructor of the `datamoments` class family
setMethod(
  "getOver", "ddatamoments",
  function(object) {
    return(object@over)
  }
)

#' Getter method of `ddatamoments` class.
#' 
#' Returns the `smoments` slot.
#' 
#' @param object An `ddatamoments` object.
#' @returns The `smoments` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' # Generate a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 0.7)), K = 2)
#' # Simulate data from the model.
#' f_data <- simulate(f_model)
#' # Calculate the mixture moments.
#' f_datamoms <- datamoments(f_data)
#' # Get the moments for the included indicators of the data. 
#' getZero(f_datamoms)
#' 
#' @seealso 
#' * [datamoments-class] for the base class for model moments
#' * [datamoments()] for the constructor of the `datamoments` class family
setMethod(
  "getZero", "ddatamoments",
  function(object) {
    return(object@zero)
  }
)

## Setters ##
## No setters as users should not manipulate a 'ddatamoments' object ##

### Private functions
### These functions are not exported
#' Generate data moments for finite mixture data
#' 
#' @description 
#' Only called implicitly. generates all moments of finite mixture data in a 
#' `fdata` object.
#' 
#' @param object A `ddatamoments` object to contain all calculated
#'   moments. 
#' @returns A `ddatamoments` object containing all moments of the 
#'   inite mixture data.
#' @noRd
".generateDdatamoments" <- function(object) {
  ## enforce column-wise ordering ##
  hasY(object@fdata, verbose = TRUE)
  datam <- getColY(object@fdata)
  ## Compute factorial moments ##
  ## fact.moments is a L x r array (L = 4) ##
  momentsm <- array(NA, dim = c(4, object@fdata@r))
  means <- apply(datam, 2, mean, na.rm = TRUE)
  object@mean <- means
  object@var <- var(datam, na.rm = TRUE)
  momentsm[1, ] <- t(means)
  momentsm[2, ] <- apply(datam * apply(datam - 1, 2, max, 0),
    2, mean,
    na.rm = TRUE
  )
  momentsm[3, ] <- apply(datam * apply(datam - 2, 2, max, 0),
    2, mean,
    na.rm = TRUE
  )
  momentsm[4, ] <- apply(datam * apply(datam - 3, 2, max, 0),
    2, mean,
    na.rm = TRUE
  )
  dimnames(momentsm) <- list(
    c("1st", "2nd", "3rd", "4th"),
    colnames(datam)
  )
  object@factorial <- momentsm
  ## Overdispersions and fractions of zeros ##
  ## over and zeros are r x 1 matrices ##
  object@over <- diag(var(datam)) - means
  object@zero <- apply(apply(datam, 2, "==", 0), 2, sum)
  return(object)
}
