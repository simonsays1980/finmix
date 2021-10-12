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

#' Finmix `mcmcpermind` class
#' 
#' @description 
#' This class defines objects to store the outputs from permuting the MCMC 
#' samples. Due to label switching the sampled component parameters are usually 
#' not assigned to the same component in each iteration. To overcome this issue 
#' the samples are permuted by using a relabeling algorithm (usually K-means) 
#' to reassign parameters. Note that due to assignment of parameters from the 
#' same iteration to the same component, the sample size could shrink. 
#' 
#' This class stores the permuted parameters together with the new MCMC sample 
#' size and the mixture log-likelihood, the prior log-likelihood, and the 
#' complete data posterior log-likelihood. All this slots are inherited from 
#' the parent class `mcmcpermfix`. In addition to these slots this class adds 
#' permuted weights, permuted indicators as well as the entropies and number 
#' of assigned observations per component.
#' 
#' @slot relabel A character defining the used algorithm for relabeling.
#' @slot weightperm An array of dimension `Mperm x K` containing the 
#'   relabeled weight parameters. 
#' @slot entropyperm An `array` of dimension `Mperm x 1` containing the 
#'   entropy for each MCMC permuted draw.
#' @slot STperm An `array` of dimension `Mperm x 1` containing all permuted 
#'   MCMC states, for the last observation in slot `@@y` of the `fdata` object 
#'   passed in to [mixturemcmc()] where a state is defined for non-Markov 
#'   models as the last indicator of this observation. 
#' @slot Sperm An `array` of dimension `N x storeS` containing the last 
#'   `storeS` permuted indicators. `storeS` is defined in the slot `@@storeS` 
#'   of the `mcmc` object passed into [mixturemcmc()].
#' @slot NKperm An `array` of dimension `Mperm x K` containing the numbers 
#'   of observations assigned to each component.
#' @exportClass mcmcpermind
#' @describeIn mcmcperm_class
#' 
#' @seealso 
#' * [mcmcpermute()] for the calling function
#' * [mcmcperfix][mcmcperm_class] for the corresponding class for models with 
#'   fixed indicators
.mcmcpermind <- setClass("mcmcpermind",
  representation(
    relabel = "character",
    weightperm = "array",
    entropyperm = "array",
    STperm = "array",
    Sperm = "array",
    NKperm = "array"
  ),
  contains = c("mcmcpermfix"),
  validity = function(object) {
    ## else: OK
    TRUE
  },
  prototype(
    relabel = character(),
    weightperm = array(),
    entropyperm = array(),
    STperm = array(),
    Sperm = array(),
    NKperm = array()
  )
)

## Getters ##

#' Getter method of `mcmcpermind` class.
#' 
#' Returns the `relabel` slot.
#' 
#' @param object An `mcmcpermind` object.
#' @returns The `relabel` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' \dontrun{getRelabel(mcmcperm)}
#' 
#' @seealso 
#' * [mcmcoutputpermbase-class] for the inheriting class
#' * [mcmcpermute()] for function permuting MCMC samples
setMethod(
  "getRelabel", "mcmcpermind",
  function(object) {
    return(object@relabel)
  }
)

#' Getter method of `mcmcpermind` class.
#' 
#' Returns the `weightperm` slot.
#' 
#' @param object An `mcmcpermind` object.
#' @returns The `weightperm` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' \dontrun{getWeightperm(mcmcperm)}
#' 
#' @seealso 
#' * [mcmcoutputpermbase-class] for the inheriting class
#' * [mcmcpermute()] for function permuting MCMC samples
setMethod(
  "getWeightperm", "mcmcpermind",
  function(object) {
    return(object@weightperm)
  }
)

#' Getter method of `mcmcpermind` class.
#' 
#' Returns the `entropyperm` slot.
#' 
#' @param object An `mcmcpermind` object.
#' @returns The `entropyperm` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' \dontrun{getEntropyperm(mcmcperm)}
#' 
#' @seealso 
#' * [mcmcoutputpermbase-class] for the inheriting class
#' * [mcmcpermute()] for function permuting MCMC samples
setMethod(
  "getEntropyperm", "mcmcpermind",
  function(object) {
    return(object@entropyperm)
  }
)

#' Getter method of `mcmcpermind` class.
#' 
#' Returns the `STperm` slot.
#' 
#' @param object An `mcmcpermind` object.
#' @returns The `STperm` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' \dontrun{getSTperm(mcmcperm)}
#' 
#' @seealso 
#' * [mcmcoutputpermbase-class] for the inheriting class
#' * [mcmcpermute()] for function permuting MCMC samples
setMethod(
  "getSTperm", "mcmcpermind",
  function(object) {
    return(object@STperm)
  }
)

#' Getter method of `mcmcpermind` class.
#' 
#' Returns the `Sperm` slot.
#' 
#' @param object An `mcmcpermind` object.
#' @returns The `Sperm` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' \dontrun{getSperm(mcmcperm)}
#' 
#' @seealso 
#' * [mcmcoutputpermbase-class] for the inheriting class
#' * [mcmcpermute()] for function permuting MCMC samples
setMethod(
  "getSperm", "mcmcpermind",
  function(object) {
    return(object@STperm)
  }
)

#' Getter method of `mcmcpermind` class.
#' 
#' Returns the `NKperm` slot.
#' 
#' @param object An `mcmcpermind` object.
#' @returns The `NKperm` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' \dontrun{getNKperm(mcmcperm)}
#' 
#' @seealso 
#' * [mcmcoutputpermbase-class] for the inheriting class
#' * [mcmcpermute()] for function permuting MCMC samples
setMethod(
  "getNKperm", "mcmcpermind",
  function(object) {
    return(object@STperm)
  }
)

## No setters as users are not intended to modify these ##
## obect.                                               ##
