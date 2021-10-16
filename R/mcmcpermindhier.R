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

#' Finmix `mcmcpermindhier` class
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
#' permuted parameter samples of the hierarchical prior.
#' 
#' Note that for models with fixed indicators `weight`s do not get permuted.
#' 
#' @slot hyperperm A named list containing the (permuted) parameters of the 
#'   hierarchical prior. 
#' @exportClass mcmcpermindhier
#' @rdname mcmcpermindhier-class
#' @keywords internal
#' 
#' @seealso 
#' * [mcmcpermute()] for the calling function
#' * [mcmcpermind-class] for the parent class definition
#' * [mcmcpermfixhier-class] for the corresponding class for models with 
#'   fixed indicators
.mcmcpermindhier <- setClass("mcmcpermindhier",
                             representation(hyperperm = "list"),
                             contains = c("mcmcpermind"),
                             validity = function(object) {
                               ## else: OK
                               TRUE
                             },
                             prototype(hyperperm = list())
)

## Getters ##

#' Getter method of `mcmcperminfhier` class.
#' 
#' Returns the `hyperperm` slot.
#' 
#' @param object An `mcmcpermindhier` object.
#' @returns The `hyperperm` slot of the `object`.
#' @exportMethod getHyperperm
#' @noRd
#' 
#' @examples 
#' \dontrun{getHyperpem(mcmcperm)}
#' 
#' @seealso 
#' * \code{\link{mcmcoutputpermind}} for the inheriting class
#' * \code{\link{mcmcpermute}} for function permuting MCMC samples
setMethod(
  "getHyperperm", "mcmcpermfixpost",
  function(object) {
    return(object@hyperperm)
  }
)
## No setters implemented as users are not intended to
## manipulate this object
