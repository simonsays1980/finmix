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

## 'datamoments' is a virtual class from which the corresponding
## datamoments for 'continuous' and 'discrete' inherit

#' Finmix `datamoments` class
#' 
#' Stores moments of a corresponding `fdata` object. 
#' 
#' @slot mean A numeric storing the mean of the slot `y` in the `fdata` object.
#' @slot var A matrix storing the variance(s and covariances) of the `y` slot 
#'   in the `fdata` object.
#' @slot VIRTUAL Virtual class containing further data moments.
#' @exportClass datamoments
#' @rdname datamoments-class
#' @seealso 
#' * [cdatamoments-class] for data moments of continuous data
#' * [ddatamoments-class] for data moments of discrete data
#' * [sdatamoments-class] for data moments of the indicators
#' 
.datamoments <- setClass(
  "datamoments",
  representation(
    mean = "numeric",
    var = "matrix",
    fdata = "fdata",
    "VIRTUAL"
  )
)


## mutual constructor for all type of datamoments ##
#' Constructor for `datamoments` classes
#' 
#' @description
#' Calling [datamoments()] generates the datamoments for an `fdata` object. 
#' Depending on the type of data either an `cdatamoments` or `ddatamoments` 
#' object is generated. If in addition the `fdata` object containes fixed 
#' indicators, these `datamoments` object also hold an `sdatamoments` class to
#' store the data moments of these indicators. 
#' 
#' @param value An `fdata` object with at least slot `y` non-empty. 
#' @returns An `datamoments` object containing the data moments for slot `y` 
#' and if available slot `S`. 
#' @export
#' 
#' @examples 
#' # Create an fdata class with Poisson data.
#' f_data <- fdata(rpois(100, 312), sim=TRUE)
#' # Compute the data moments.
#' datamoments(f_data)
#' 
#' @seealso 
#' * [datamoments-class] for all slots of this class
#' * [cdatamoments-class] for the class for continuous data
#' * [ddatamoments-class] for the class for discrete data
"datamoments" <- function(value = fdata()) {
  hasY(value, verbose = TRUE)
  if (value@type == "continuous") {
    .Object <- .cdatamoments(value = value)
  } else {
    .Object <- .ddatamoments(value = value)
  }
  return(.Object)
}
