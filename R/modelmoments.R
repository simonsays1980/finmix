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
#' Finmix `modelmoments` class
#' 
#' @description
#' Defines a container to hold the moments of a finite mixture model. The 
#' finmix `model` object should contains parameters and weights. 
#' 
#' @slot mean A vector of component means. 
#' @slot var An array of components variances or in case of multivariate 
#' distributions covariance matrices. 
#' @slot model The corresponding `model` object.
#' @exportClass modelmoments
#' 
#' @name modelmoments_class
#' @seealso 
#' * [modelmoments()] the constructor of the `modelmoments` class
setClass("modelmoments",
  representation(
    mean = "vector",
    var = "array",
    model = "model"
  ),
  validity = function(object) {
    ## else: OK
    TRUE
  },
  prototype(
    mean     = vector(),
    var      = array(),
    model    = model()
  )
)

#' Constructor of finmix `modelmoments` class
#' 
#' Calling [modelmoments()] calculates the corresponding moments of the 
#' finite mixture distribution defined in the `model` object. The `model` 
#' object should contain parameters in slot `par` and weights in slot `weight`.
#' 
#' @param model A `model` object containing defined parameters in slot `par` 
#'   and defined weights in slot `weight`.
#' @returns A `modelmoments` object with calculated moments of the finite 
#' mixture model defined in the `model` object. 
#' @export
#' 
#' @examples 
#' f_model <- model("poisson", par=list(lambda=c(0.3, 0.1)), 
#'                  weight=matrix(c(0.3, 0.7), nrow=1))
#' modelmoments(f_model)
#' 
#' @seealso 
#' * [modelmoments_class] for all slots of the `modelmoments` class
"modelmoments" <- function(model) {
  dist <- model@dist
  if (dist == "normult") {
    .normultmodelmoments(model = model)
  } else if (dist == "studmult") {
    .studmultmodelmoments(model = model)
  } else if (dist == "student") {
    .studentmodelmoments(model = model)
  } else if (dist == "normal") {
    .normalmodelmoments(model = model)
  } else if (dist == "exponential") {
    .exponentialmodelmoments(model = model)
  } else if (dist %in% c("poisson", "cond.poisson")) {
    .poissonmodelmoments(model = model)
  } else if (dist == "binomial") {
    .binomialmodelmoments(model = model)
  }
}

## Getters ##
#' Getter method of `modelmoments` class.
#' 
#' Returns the `mean` slot of a `modelmoments` object. 
#' 
#' @param object A `modelmoments` object.
#' @returns The `mean` slot of the `object`.
#' @exportMethod getMean
#' @describeIn modelmoments_class
#' 
#' @examples 
#' f_model <- model("poisson", par=list(lambda=c(0.3, 0.1)), 
#'                  weight=matrix(c(0.3, 0.7), nrow=1))
#' f_moments <- modelmoments(f_model)
#' getMean(f_moments)
#' 
#' @seealso [modelmoments_class] for all slots of the `modelmoments` class
setMethod(
  "getMean", "modelmoments",
  function(object) {
    return(object@mean)
  }
)

#' Getter method of `modelmoments` class.
#' 
#' Returns the `var` slot of a `modelmoments` object. 
#' 
#' @param object A `modelmoments` object.
#' @returns The `var` slot of the `object`.
#' @exportMethod getVar
#' @describeIn modelmoments_class
#' 
#' @examples 
#' f_model <- model("poisson", par=list(lambda=c(0.3, 0.1)), 
#'                  weight=matrix(c(0.3, 0.7), nrow=1))
#' f_moments <- modelmoments(f_model)
#' getVar(f_moments)
#' 
#' @seealso [modelmoments_class] for all slots of the `modelmoments` class
setMethod(
  "getVar", "modelmoments",
  function(object) {
    return(object@var)
  }
)

#' Getter method of `modelmoments` class.
#' 
#' Returns the `model` slot of a `modelmoments` object. 
#' 
#' @param object A `modelmoments` object.
#' @returns The `model` slot of the `object`.
#' @exportMethod getModel
#' @describeIn modelmoments_class
#' 
#' @examples 
#' f_model <- model("poisson", par=list(lambda=c(0.3, 0.1)), 
#'                  weight=matrix(c(0.3, 0.7), nrow=1))
#' f_moments <- modelmoments(f_model)
#' getModel(f_moments)
#' 
#' @seealso [modelmoments_class] for all slots of the `modelmoments` class
setMethod(
  "getModel", "modelmoments",
  function(object) {
    return(object@model)
  }
)

## Setters are not provided as users are not intended to manipulate ##
## this object ##
