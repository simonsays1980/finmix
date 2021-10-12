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

#' Finmix `mcmcpermfix` class
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
#' complete data posterior log-likelihood.
#' 
#' Note that for models with fixed indicators `weight`s do not get permuted.
#' 
#' @slot Mperm An integer storing the MCMC sample size after relabeling.
#' @slot parperm A named list containing the permuted component parameters. 
#' @slot logperm A named list containing the mixture log-likelihood, the prior 
#'   log-likelihood, and the complete data posterior log-likelihood.
#' @exportClass mcmcpermfix
#' @family mcmcperm-classes
#' @rdname mcmcpermfix-class
#' 
#' @seealso 
#' * [mcmcpermute()] for the calling function
#' * [mcmcpermind-class] for the corresponding class for models with 
#'   unknown indicators
.mcmcpermfix <- setClass("mcmcpermfix",
  representation(
    Mperm       = "integer",
    parperm     = "list",
    logperm     = "list"
  ),
  validity = function(object) {
    ## else: OK
    TRUE
  },
  prototype(
    Mperm = integer(),
    parperm = list(),
    logperm = list()
  )
)

## Getters ##

#' Getter method of `mcmcpermfix` class.
#' 
#' Returns the `Mperm` slot.
#' 
#' @param object An `mcmcpermfix` object.
#' @returns The `Mperm` slot of the `object`.
#' @aliases mcmcpermfix_class, mcmcpermfixhier_class, mcmcpermfixpost_class, 
#'   mcmcpermfixhierpost
#' 
#' @examples 
#' \dontrun{getMperm(mcmcperm)}
#' 
#' @seealso 
#' * [mcmcoutputpermfix-class] for the inheriting class
#' * [mcmcpermute()] for function permuting MCMC samples
setMethod(
  "getMperm", "mcmcpermfix",
  function(object) {
    return(object@Mperm)
  }
)

#' Getter method of `mcmcpermfix` class.
#' 
#' Returns the `parperm` slot.
#' 
#' @param object An `mcmcpermfix` object.
#' @returns The `parperm` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' \dontrun{getParperm(mcmcperm)}
#' 
#' @seealso 
#' * [mcmcoutput-class] for the inheriting class
#' * [mcmcpermute()] for function permuting MCMC samples
setMethod(
  "getParperm", "mcmcpermfix",
  function(object) {
    return(object@parperm)
  }
)

#' Getter method of `mcmcpermfix` class.
#' 
#' Returns the `logperm` slot.
#' 
#' @param object An `mcmcpermfix` object.
#' @returns The `logperm` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' \dontrun{getLogperm(mcmcperm)}
#' 
#' @seealso 
#' * [mcmcoutput-class] for the inheriting class
#' * [mcmcpermute()] for function permuting MCMC samples
setMethod(
  "getLogperm", "mcmcpermfix",
  function(object) {
    return(object@logperm)
  }
)

## No setters as users are not intended to modify these ##
## objects.                                             ##
