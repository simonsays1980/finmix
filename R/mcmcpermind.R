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
setMethod(
  "getRelabel", "mcmcpermind",
  function(object) {
    return(object@relabel)
  }
)

setMethod(
  "getWeightperm", "mcmcpermind",
  function(object) {
    return(object@weightperm)
  }
)

setMethod(
  "getEntropyperm", "mcmcpermind",
  function(object) {
    return(object@entropyperm)
  }
)

setMethod(
  "getSTperm", "mcmcpermind",
  function(object) {
    return(object@STperm)
  }
)

setMethod(
  "getSperm", "mcmcpermind",
  function(object) {
    return(object@STperm)
  }
)

setMethod(
  "getNKperm", "mcmcpermind",
  function(object) {
    return(object@STperm)
  }
)

## No setters as users are not intended to modify these ##
## obect.                                               ##
