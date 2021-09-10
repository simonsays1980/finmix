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

.dmodelmoments <- setClass("dmodelmoments",
  representation(
    over            = "numeric",
    factorial       = "array",
    zero            = "numeric"
  ),
  contains = c("modelmoments"),
  validity = function(object) {
    ## else: OK
    TRUE
  },
  prototype(
    over       = numeric(),
    factorial  = array(),
    zero       = numeric()
  )
)

## Getters ##
setMethod("getOver", "dmodelmoments", function(object) {
  return(object@over)
})

setMethod("getFactorial", "dmodelmoments", function(object) {
  return(object@factorial)
})

setMethod("getZero", "dmodelmoments", function(object) {
  return(object@zero)
})

## Setters ##
## No setters as users should not manipulate a 'dmodelmoments' object ##
