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

#' Finmix `cmodelmoments` class
#' 
#' @description 
#' This class defines the general theoretical moments of a finite mixture model 
#' with continuous data. 
#' 
#' @slot higher An array containing the four higher centralized moments of the 
#'   (in case of multivariate data marginal) finite mixture.
#' @slot skewness A vector containing the skewness(es) of the finite mixture 
#'   model.
#' @slot kurtosis A vector containing the kurtosis(es) of the finite mixture 
#'   model.
#' @exportClass cmodelmoments
#' @name cmodelmoments
#' 
#' @seealso 
#' * [modelmoments] for the base class
#' * [modelmoments()] for the constructor of any `modelmoments` inherited class
.cmodelmoments <- setClass("cmodelmoments",
  representation(
    higher      = "array",
    skewness    = "vector",
    kurtosis    = "vector"
  ),
  contains = c("modelmoments"),
  validity = function(object) {
    ## else: OK
    TRUE
  },
  prototype(
    higher     = array(),
    skewness   = vector(),
    kurtosis   = vector()
  )
)

## Getters ##
#' Getter method of `cmodelmoments` class.
#' 
#' Returns the `higher` slot.
#' 
#' @param object An `cmodelmoments` object.
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
setMethod("getHigher", "cmodelmoments", function(object) {
  return(object@higher)
})

#' Getter method of `cmodelmoments` class.
#' 
#' Returns the `skewness` slot.
#' 
#' @param object An `cmodelmoments` object.
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
setMethod("getSkewness", "cmodelmoments", function(object) {
  return(object@skewness)
})

#' Getter method of `cmodelmoments` class.
#' 
#' Returns the `kurtosis` slot.
#' 
#' @param object An `cmodelmoments` object.
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
setMethod("getKurtosis", "cmodelmoments", function(object) {
  return(object@kurtosis)
})

## Setters ##
## No setters as users should not manipulate a 'nsmodelmoments' object ##
