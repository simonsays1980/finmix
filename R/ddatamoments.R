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

.ddatamoments <- setClass("ddatamoments",
                          representation(factorial   = "array",
                                         over        = "vector",
                                         zero        = "vector",
                                         smoments    = "sdatamomentsOrNULL"),
                          contains = c("datamoments"),
                          validity = function(object) 
                          {
                              ## else: ok
                              TRUE
                          },
                          prototype(factorial   = array(),
                                    over        = vector(),
                                    zero        = vector(),
                                    smoments    = .sdatamoments()
                                    )
)

setMethod("initialize", "ddatamoments", 
          function(.Object, ..., value = fdata()) 
          {
              .Object@fdata <- value
			  if (hasS(value)) {
                  .Object@smoments  <- sdatamoments(value)
			  } else {
                  .Object@smoments  <- NULL
              }
              generateMoments(.Object)
          }
)

## Generic set in 'groupmoments.R' ##
setMethod("generateMoments", "ddatamoments",
          function(object) 
          {
              .generateDdatamoments(object)
          }
)

setMethod("show", "ddatamoments", 
          function(object) 
          {
              cat("Object 'datamoments'\n")
              cat("     mean        : Vector of", 
                  length(object@mean), "\n")
              cat("     var         : Vector of",
                  length(object@var), "\n")
              cat("     factorial   :",
                  paste(dim(object@factorial), collapse = "x"), "\n")
              cat("     over        : Vector of",
                  length(object@over), "\n")
              cat("     zero        : Vector of",
                  length(object@zero), "\n")
              if (hasS(object@fdata)) {
                  cat("     smoments    : Object of class",
                      class(object@smoments), "\n")
              }
              cat("     fdata       : Object of class",
                  class(object@fdata), "\n")
          }
)

## Getters ##
setMethod("getSmoments", "ddatamoments", 
          function(object) 
          {
              return(object@smoments)
          }
)

setMethod("getFactorial", "ddatamoments", 
          function(object) 
          {
              return(object@factorial)
          }
)

setMethod("getOver", "ddatamoments", 
          function(object) 
          {
              return(object@over)
          }
)

setMethod("getZero", "ddatamoments", 
          function(object) 
          {
              return(object@zero)
          }
)

## Setters ##
## No setters as users should not manipulate a 'ddatamoments' object ##

### Private functions
### These functions are not exported
".generateDdatamoments" <- function(object) 
{
    ## enforce column-wise ordering ##
    hasY(object@fdata, verbose = TRUE)
    datam   <- getColY(object@fdata)
    ## Compute factorial moments ##
    ## fact.moments is a L x r array (L = 4) ## 
    momentsm <- array(NA, dim = c(4, object@fdata@r))
    means           <- apply(datam, 2, mean, na.rm = TRUE)
    object@mean <- means
    object@var    <- var(datam, na.rm = TRUE)
    momentsm[1, ]   <- t(means)  
    momentsm[2, ]   <- apply(datam * apply(datam - 1, 2, max, 0), 
                             2, mean, na.rm = TRUE)
    momentsm[3, ]   <- apply(datam * apply(datam - 2, 2, max, 0),
                             2, mean, na.rm = TRUE)
    momentsm[4, ]   <- apply(datam * apply(datam - 3, 2, max, 0),
                             2, mean, na.rm = TRUE)
    dimnames(momentsm) <- list(c("1st", "2nd", "3rd", "4th"),
                               colnames(datam))
    object@factorial <- momentsm
    ## Overdispersions and fractions of zeros ##
    ## over and zeros are r x 1 matrices ## 
    object@over <- diag(var(datam)) - means
    object@zero <- apply(apply(datam, 2, "==", 0), 2, sum)
    return(object)
}
