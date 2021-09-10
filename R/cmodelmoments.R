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
setMethod("getHigher", "cmodelmoments", function(object) {
  return(object@higher)
})

setMethod("getSkewness", "cmodelmoments", function(object) {
  return(object@skewness)
})

setMethod("getKurtosis", "cmodelmoments", function(object) {
  return(object@kurtosis)
})

## Setters ##
## No setters as users should not manipulate a 'nsmodelmoments' object ##
