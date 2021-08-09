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

.exponentialmodelmoments <- setClass("exponentialmodelmoments",
                                     representation(B   = "numeric",
                                                    W   = "numeric",
                                                    R   = "numeric"),
                                     contains = c("cmodelmoments"),
                                     validity = function(object) {
                                         ## else: OK
                                         TRUE
                                     },
                                     prototype(B    = numeric(),
                                               W    = numeric(),
                                               R    = numeric()
                                               )
)

setMethod("initialize", "exponentialmodelmoments",
          function(.Object, ..., model) 
          {
              .Object <- callNextMethod(.Object, ..., model = model)
              generateMoments(.Object)
          }
)

setMethod("generateMoments", "exponentialmodelmoments",
          function(object) 
          {
              .generateMomentsExponential(object)
          }
)

setMethod("show", "exponentialmodelmoments",
          function(object) {
              cat("Object 'modelmoments'\n")
              cat("     mean        : Vector of", 
                  length(object@mean), "\n")
              cat("     var         :", 
                  paste(dim(object@var), collapse = "x"), "\n")
              cat("     higher      :", 
                  paste(dim(object@higher), collapse = "x"),
                  "\n")
              cat("     skewness    :", object@skewness, "\n")
              cat("     kurtosis    :", object@kurtosis, "\n")
              cat("     B           :", object@B, "\n")
              cat("     W           :", object@W, "\n")
              cat("     R           :", object@R, "\n")
              cat("     model       : Object of class", 
                  class(object@model), "\n")
          }
)

setMethod("getB", "exponentialmodelmoments",
          function(object) 
          {
              return(object@B)
          }
)

setMethod("getW", "exponentialmodelmoments",
          function(object) 
          {
              return(object@W)
          }
)

setMethod("getR", "exponentialmodelmoments",
          function(object)
          {
              return(object@R)
          }
)

## No setters as users are not intended to manipulate ##
## this object ##

### Private functions 
### These functions are not exported 
".generateMomentsExponential" <- function(object) 
{
    lambda          <- object@model@par$lambda
    weight          <- object@model@weight
    object@mean     <- sum(weight * 1/lambda)
    highm <- .mixturemoments.exponential(object@model, 4, object@mean)
    dimnames(highm) <- list(c("1st", "2nd", "3rd", "4th"), "")
    object@higher   <- highm
    object@var      <- array(object@higher[2], dim = c(1, 1))
    object@skewness <- object@higher[3]/object@higher[2]^1.5
    object@kurtosis <- object@higher[4]/object@higher[2]^2
    object@W        <- sum(weight * 1/lambda^2)
    object@B        <- sum(weight * (1/lambda - object@mean)^2)
    object@R        <- 1 - object@W/object@var[1]
    return(object)
}

