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
# along with Rcpp.  If not, see <http://www.gnu.org/licenses/>.

#' Finmix `mcmcpermfixpost` class
#' 
#' @description 
#' This class defines objects to store the outputs from permuting the MCMC 
#' samples. Due to label switching the sampled component parameters are usually 
#' not assigned to the same component in each iteration. To overcome this issue 
#' the samples are permuted by using a relabeling algorithm (usually K-means) 
#' to reassign parameters. Note that due to assignment of parameters from the 
#' same iteration to the same component, the sample size could shrink. 
#' 
#' This class is supplementing the parent class by adding a slot to store the 
#' permuted parameter samples of the posterior densities.
#' 
#' Note that for models with fixed indicators `weight`s do not get permuted.
#' 
#' @slot postperm A named list containing a named list `par` with array(s) of 
#'   parameters from the posterior density. 
#'   
#' @exportClass mcmcpermfixpost
#' @rdname mcmcpermfixpost-class
#' 
#' @seealso 
#' * [mcmcpermute()] for the calling function
#' * [mcmcpermfix-class] for the parent class definition
#' * [mcmcpermindpost-class]for the corresponding class for models with 
#'   unknown indicators
.mcmcpermfixpost <- setClass("mcmcpermfixpost",
  representation(postperm = "list"),
  contains = c("mcmcpermfix"),
  validity = function(object) {
    ## else: OK
    TRUE
  },
  prototype(postperm = list())
)

## Getters ##

#' Getter method of `mcmcpermfixpost` class.
#' 
#' Returns the `postperm` slot.
#' 
#' @param object An `mcmcpermfixpost` object.
#' @returns The `postperm` slot of the `object`.
#' @exportMethod getPostperm
#' @keywords internal
#' 
#' @examples 
#' \dontrun{getMperm(mcmcperm)}
#' 
#' @seealso 
#' * [mcmcoutputpermfix-class] for the inheriting class
#' * [mcmcpermute()] for function permuting MCMC samples
setMethod(
  "getPostperm", "mcmcpermfixpost",
  function(object) {
    return(object@postperm)
  }
)

## No setters implemented as users are not intended to
## manipulate this object
