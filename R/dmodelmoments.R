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
#' Finmix `dmodelmoments` class
#' 
#' @description 
#' This class defines the general theoretical moments of a finite mixture model 
#' with discrete data. 
#' 
#' @slot over A numeric containing the over-dispersion.
#' @slot factorial An array containing the first four factorial moments.
#' @slot zero An numeric cotaining the excess zeros.  
#' @exportClass dmodelmoments
#' @name dmodelmoments
#' 
#' @seealso 
#' * [modelmoments] for the base class
#' * [modelmoments()] for the constructor of any `modelmoments` inherited class
.dmodelmoments <- setClass("dmodelmoments",
  representation(
    over            = "numeric",
    factorial       = "array",
    zero            = "numeric"
  ),
  contains = c("modelmoments"),
  validity = function(object) {
    ## else: OK
    TRUE
  },
  prototype(
    over       = numeric(),
    factorial  = array(),
    zero       = numeric()
  )
)

## Getters ##
#' Getter method of `dmodelmoments` class.
#' 
#' Returns the `higher` slot.
#' 
#' @param object An `dmodelmoments` object.
#' @returns The `higher` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' f_model <- model("c", par=list(lambda=c(0.3, 0.1)), 
#'                  weight=matrix(c(0.3, 0.7), nrow=1))
#' f_moments <- modelmoments(f_model)
#' getHigher(f_moments)
#' 
#' @seealso 
#' * [modelmoments] for the base class for model moments
#' * [modelmoments()] for the constructor of the `modelmoments` class family
setMethod("getOver", "dmodelmoments", function(object) {
  return(object@over)
})

#' Getter method of `dmodelmoments` class.
#' 
#' Returns the `skewness` slot.
#' 
#' @param object An `dmodelmoments` object.
#' @returns The `skewness` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' f_model <- model("c", par=list(lambda=c(0.3, 0.1)), 
#'                  weight=matrix(c(0.3, 0.7), nrow=1))
#' f_moments <- modelmoments(f_model)
#' getSkewness(f_moments)
#' 
#' @seealso 
#' * [modelmoments] for the base class for model moments
#' * [modelmoments()] for the constructor of the `modelmoments` class family
setMethod("getFactorial", "dmodelmoments", function(object) {
  return(object@factorial)
})

#' Getter method of `dmodelmoments` class.
#' 
#' Returns the `kurtosis` slot.
#' 
#' @param object An `dmodelmoments` object.
#' @returns The `kurtosis` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' f_model <- model("c", par=list(lambda=c(0.3, 0.1)), 
#'                  weight=matrix(c(0.3, 0.7), nrow=1))
#' f_moments <- modelmoments(f_model)
#' getKurtosis(f_moments)
#' 
#' @seealso 
#' * [modelmoments] for the base class for model moments
#' * [modelmoments()] for the constructor of the `modelmoments` class family
setMethod("getZero", "dmodelmoments", function(object) {
  return(object@zero)
})

## Setters ##
## No setters as users should not manipulate a 'dmodelmoments' object ##
