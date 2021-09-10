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

.sdatamoments <- setClass("sdatamoments",
  representation(
    gmoments = "groupmoments",
    fdata = "fdata"
  ),
  validity = function(object) {
    ## else: OK
    TRUE
  }
)

setClassUnion("sdatamomentsOrNULL", members = c("sdatamoments", "NULL"))

## mutual constructor for both types of sdatamoments ##
"sdatamoments" <- function(value = fdata()) {
  hasY(value, verbose = TRUE)
  hasS(value, verbose = TRUE)
  if (value@type == "discrete") {
    object <- .sdatamoments(value = value)
  } else {
    object <- .csdatamoments(value = value)
  }
  return(object)
}

setMethod(
  "initialize", "sdatamoments",
  function(.Object, ..., value = fdata()) {
    .Object@fdata <- value
    .Object@gmoments <- .groupmoments(value = value)
    return(.Object)
  }
)

setMethod(
  "show", "sdatamoments",
  function(object) {
    cat("Object 'sdatamoments'\n")
    cat(
      "     gmoments    : Object of class",
      class(object@gmoments), "\n"
    )
    cat(
      "     fdata       : Object of class",
      class(object@fdata), "\n"
    )
  }
)

## Getters ##
setMethod(
  "getGmoments", "sdatamoments",
  function(object) {
    return(object@gmoments)
  }
)

setMethod(
  "getFdata", "sdatamoments",
  function(object) {
    return(object@fdata)
  }
)

## Setters ##
## No Setters, as it is adviced for users not to manipulate moment objects ##
