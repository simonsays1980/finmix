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

.binomialmodelmoments <- setClass("binomialmodelmoments",
  representation(extrabinvar = "numeric"),
  contains = c("dmodelmoments"),
  validity = function(object) {
    ## else: OK
    TRUE
  },
  prototype(extrabinvar = numeric())
)

setMethod(
  "initialize", "binomialmodelmoments",
  function(.Object, ..., model) {
    .Object <- callNextMethod(.Object, ..., model = model)
    generateMoments(.Object)
  }
)

setMethod(
  "generateMoments", "binomialmodelmoments",
  function(object) {
    .generateMomentsBinomial(object)
  }
)

setMethod(
  "show", "binomialmodelmoments",
  function(object) {
    cat("Object 'modelmoments'\n")
    cat(
      "     mean        : Vector of",
      length(object@mean), "\n"
    )
    cat(
      "     var         :",
      paste(dim(object@var), collapse = "x"), "\n"
    )
    cat(
      "     factorial   :",
      paste(dim(object@factorial), collapse = "x"),
      "\n"
    )
    cat("     over        :", object@over, "\n")
    cat("     zero        :", object@zero, "\n")
    cat(
      "     extrabinvar :", object@extrabinvar,
      "\n"
    )
    cat(
      "     model       : Object of class",
      class(object@model), "\n"
    )
  }
)

## Getters ##
setMethod(
  "getExtrabinvar", "binomialmodelmoments",
  function(object) {
    return(object@extrabinvar)
  }
)

## No setters as users are not intended to manipulate ##
## this object ##

### Private functions
### These function are not exported
".generateMomentsBinomial" <- function(object) {
  p <- object@model@par$p
  n <- object@model@par$n
  weight <- object@model@weight
  object@mean <- sum(weight * n * p)
  object@var <- array(sum(weight * (n * p - object@mean)^2)
  + sum(weight * n * p * (1 - p)), dim = c(1, 1))
  factm <- array(NA, dim = c(4, 1))
  factm[1] <- object@mean
  for (i in seq(2, 4)) {
    if (n >= i) {
      factm[i] <- sum(weight * factorial(n) / factorial(n - i) * p^i)
    } else {
      factm[i] <- NaN
    }
  }
  dimnames(factm) <- list(c("1st", "2nd", "3rd", "4th"), "")
  object@factorial <- factm
  if (object@model@K > 1) {
    object@over <- object@var[1] - object@mean
  } else {
    object@over <- 0
  }
  object@zero <- sum(weight * (1 - p)^n)
  object@extrabinvar <- object@mean * (1 - object@mean / n[1])
  return(object)
}
