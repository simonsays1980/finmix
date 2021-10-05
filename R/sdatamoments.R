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

#' Finmix `sdatamoments` class
#' 
#' Stores moments for indicators of discrete data. 
#' 
#' @slot gmoments A [groupmoments][groupmoments_class] object storing the 
#'   moments for each mixture component. 
#' @slot fdata An [fdata][fdata_class] object with data from a discrete valued
#'   mixture distribution.
#' @exportClass sdatamoments
#' @name sdatamoments_class
#' @seealso 
#' * [datamoments][datamoments_class] for the base class for data moments
#' * [datamoments()] for the constructor of any object of the `datamoments` 
#'   class family
#' * [groupmoments][groupmoments_class] for the parent class 
#' * [csdatamoments][csdatamoments_class] for the corresponding class defining
#'   moments for data from a continuous-valued finite mixture
.sdatamoments <- setClass("sdatamoments",
  representation(
    gmoments = "groupmoments",
    fdata = "fdata"
  ),
  validity = function(object) {
    ## else: OK
    TRUE
  }
)

#' Finmix class union of `sdatamoments` and `NULL`
#' 
#' @description
#' Defines a class union such that the object held by a child class can also
#' be `NULL`.
#' 
#' @export
#' @keywords internal
setClassUnion("sdatamomentsOrNULL", members = c("sdatamoments", "NULL"))

## mutual constructor for both types of sdatamoments ##
#' Finmix `sdatamoments` constructor
#' 
#' Calling [sdatamoments()] constructs an object of class `sdatamoments` or
#' `csdatamoments` depending on the `type` slot of the argument `value`. If 
#' this slot is `"discrete"` an `sdatamoments` object is returned and if the 
#' slot is `"continuous"`, a `csdatamoments` object is returned. 
#' 
#' @param value An [fdata][fdata_class] object containing the indicators for 
#'   which moments should be calculated.
#' @return If slot `type` of the argument `value` is `"discrete"` an 
#'   `sdatamoments` object is returned and if the slot is `"continuous"`, 
#'   a `csdatamoments` object is returned.
#' @export
#' @name sdatamoments
#' 
#' @example
#' # Define a model of exponential mixtures.
#' f_model <- model("exponential", par = list(lambda = c(0.3, 0.7)), K = 2)
#' # Simulate data from the model.
#' f_data <- simulate(f_model)
#' # Compute data moments for the indicators.
#' sdatamoments(f_data)
#' 
#' @seealso 
#' * [sdatamoments][sdatamoments_class] for the class of indicator 
#'   moments for discrete data
#' * [csdatamoments][csdatamoments_class] for the class of indicator moments
#'   for continuous
#' * [groupmoments][groupmoments_class] for the parent class## Copyright (C) 2013 Lars Simon Zehnder
"sdatamoments" <- function(value = fdata()) {
  hasY(value, verbose = TRUE)
  hasS(value, verbose = TRUE)
  if (value@type == "discrete") {
    object <- .sdatamoments(value = value)
  } else {
    object <- .csdatamoments(value = value)
  }
  return(object)
}

#' Initializer of the `sdatamoments` class
#' 
#' @description
#' Only used implicitly. The initializer calls the constructor for a 
#' [groupmoments][groupmoments_class] object. to generate in the initialization 
#' step the moments for a passed-in `fdata` object.
#' 
#' @param .Object An object: see the "initialize Methods" section in 
#'   [initialize].
#' @param ... Arguments to specify properties of the new object, to be passed 
#'   to `initialize()`.
#' @param model A finmix [fdata][fdata_class] object containing the observations.
#' @keywords internal
#' 
#' @seealso 
#' * [Classes_Details] for details of class definitions, and 
#' * [setOldClass] for the relation to S3 classes
setMethod(
  "initialize", "sdatamoments",
  function(.Object, ..., value = fdata()) {
    .Object@fdata <- value
    .Object@gmoments <- .groupmoments(value = value)
    return(.Object)
  }
)

#' Shows a summary of an `sdatamoments` object.
#' 
#' Calling [show()] on an `sdatamoments` object gives an overview 
#' of the moments of a finite mixture with discrete data.
#' 
#' @param object An `sdatamoments` object.
#' @returns A console output listing the slots and summary information about
#'   each of them. 
#' @exportMethod show
#' @describeIn sdatamoments_class Shows a summary of an object
setMethod(
  "show", "sdatamoments",
  function(object) {
    cat("Object 'sdatamoments'\n")
    cat(
      "     gmoments    : Object of class",
      class(object@gmoments), "\n"
    )
    cat(
      "     fdata       : Object of class",
      class(object@fdata), "\n"
    )
  }
)

## Getters ##
#' Getter method of `sdatamoments` class.
#' 
#' Returns the `gmoments` slot.
#' 
#' @param object An `sdatamoments` object.
#' @returns The `gmoments` slot of the `object`.
#' @noRd
#' 
#' @exportMethod getGmoments
#' @examples 
#' # Generate a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 0.7)), K = 2)
#' # Simulate data from the model.
#' f_data <- simulate(f_model)
#' # Calculate the mixture moments.
#' f_sdatamoms <- sdatamoments(f_data)
#' # Get the moments for the included indicators of the data. 
#' getGmoments(f_sdatamoms)
#' 
#' @seealso 
#' * [datamoments][datamoments_class] for the base class for model moments
#' * [datamoments()][datamoments] for the constructor of the `datamoments` 
#'   class family
#' * [sdatamoments][sdatamoments_class] for the class definition
#' * [sdatamoments()][sdatamoments] for the constructor of the class
setMethod(
  "getGmoments", "sdatamoments",
  function(object) {
    return(object@gmoments)
  }
)

#' Getter method of `sdatamoments` class.
#' 
#' Returns the `fdata` slot.
#' 
#' @param object An `sdatamoments` object.
#' @returns The `fdata` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' # Generate a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 0.7)), K = 2)
#' # Simulate data from the model.
#' f_data <- simulate(f_model)
#' # Calculate the mixture moments.
#' f_sdatamoms <- sdatamoments(f_data)
#' # Get the moments for the included indicators of the data. 
#' getFdata(f_sdatamoms)
#' 
#' @seealso 
#' * [datamoments][datamoments_class] for the base class for model moments
#' * [datamoments()][datamoments] for the constructor of the `datamoments` 
#'   class family
#' * [sdatamoments][sdatamoments_class] for the class definition
#' * [sdatamoments()][sdatamoments] for the constructor of the class
setMethod(
  "getFdata", "sdatamoments",
  function(object) {
    return(object@fdata)
  }
)

## Setters ##
## No Setters, as it is adviced for users not to manipulate moment objects ##