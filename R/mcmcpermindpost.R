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

#' Finmix `mcmcpermindpost` class
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
#' @slot postperm A named list containing a named list `par` with array(s) of 
#'   parameters from the posterior density. 
#' @exportClass mcmcpermindpost
#' @describeIn mcmcperm_class
#' 
#' @seealso 
#' * [mcmcpermute()] for the calling function
#' * [mcmcpermind][mcmcperm_class] for the parent class definition
#' * [mcmcpermfixpost][mcmcperm_class] for the corresponding class for models 
#'   with fixed indicators
.mcmcpermindpost <- setClass("mcmcpermindpost",
  representation(postperm = "list"),
  contains = c("mcmcpermind"),
  validity = function(object) {
    ## else: OK
    TRUE
  },
  prototype(postperm = list())
)

## Getters ##

#' Getter method of `mcmcpermindpost` class.
#' 
#' Returns the `postperm` slot.
#' 
#' @param object An `mcmcpermindpost` object.
#' @returns The `postperm` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' \dontrun{getMperm(mcmcperm)}
#' 
#' @seealso 
#' * [mcmcoutputpermpost-class] for the inheriting class
#' * [mcmcoutputpermhierpost-class] for the inheriting class with 
#'   hierarchical prior
#' * [mcmcpermute()] for function permuting MCMC samples
setMethod(
  "getPostperm", "mcmcpermindpost",
  function(object) {
    return(object@postperm)
  }
)

## No setters as users are not intended to manipulate
## this objects
