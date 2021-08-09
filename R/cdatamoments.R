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

.cdatamoments <- setClass("cdatamoments",
                          representation(higher      = "array",
                                         skewness    = "vector",
                                         kurtosis    = "vector",
                                         corr        = "matrix",
                                         smoments    = "csdatamomentsOrNULL"
                                         ),
                          contains = c("datamoments"),
                          validity = function(object) 
                          {
                              ## else: ok
                              TRUE
                          },
                          prototype(higher    = array(),
                                    skewness  = vector(),
                                    kurtosis  = vector(),
                                    corr      = matrix(),
                                    smoments  = .csdatamoments()
                                    )
)

setMethod("initialize", "cdatamoments", 
          function(.Object, ..., value = fdata()) 
          {
              .Object@fdata <- value
			  if (hasS(value)) {
                  .Object@smoments  <- sdatamoments(value = value)
              } else {
                  .Object@smoments  <- NULL
              }
              generateMoments(.Object)
          }
)

setMethod("generateMoments", "cdatamoments",
          function(object) 
          {
              .generateCdatamoments(object)
          }
)

setMethod("show", "cdatamoments", 
          function(object) 
          {
              cat("Object 'datamoments'\n")
              cat("     mean        : Vector of", 
                  length(object@mean), "\n")
              cat("     var         : Vector of", 
                  length(object@var), "\n")
              cat("     higher      :", 
                  paste(dim(object@higher), collapse = "x"), "\n")
              cat("     skewness    : Vector of", 
                  length(object@skewness), "\n")
              cat("     kurtosis    : Vector of", 
                  length(object@kurtosis), "\n")
              if (!all(is.na(object@corr))) {
                  cat("     corr        :",
                      paste(dim(object@corr), collapse = "x"), "\n")
              }
              if (hasS(object@fdata)) {
                  cat("     smoments    : Object of class", 
                      class(object@smoments), "\n")
              }
              cat("     fdata        : Object of class",
                  class(object@fdata), "\n")
          }
)

## Getters ##
setMethod("getSmoments", "cdatamoments", 
          function(object) 
          {
              return(object@smoments)
          }
)

setMethod("getHigher", "cdatamoments", 
          function(object) 
          {
              return(object@higher)
          }
)

setMethod("getSkewness", "cdatamoments", 
          function(object) 
          {
              return(object@skewness)
          }
)

setMethod("getKurtosis", "cdatamoments", 
          function(object) 
          {
              return(object@kurtosis)
          }
)

setMethod("getCorr", "cdatamoments", 
          function(object) 
          {
              return(object@corr)
          }
)
## Setters ##
## No setters as users should not manipulate a 'cdatamoments' object ##

### Private functions
### These function are not exported
".generateCdatamoments" <- function(object) 
{
    ## enforce column-wise ordering ##
    hasY(object@fdata, verbose = TRUE)
    datam   <- getColY(object@fdata)
    ## Compute higher moments ##
    ## higher.moments is a r x L matrix (L = 4) ##
    means <- apply(datam, 2, mean, na.rm = TRUE)
    object@mean <- means
    object@var <- var(datam, na.rm = TRUE)
    d <- datam - rep(means, each = nrow(datam)) 
    momentsm <- array(0, dim = c(4, object@fdata@r))
    momentsm[2,] <- apply(d^2, 2, mean, na.rm = TRUE)
    momentsm[3,] <- apply(d^3, 2, mean, na.rm = TRUE)
    momentsm[4,] <- apply(d^4, 2, mean, na.rm = TRUE)
    dimnames(momentsm) <- list(c("1st", "2nd", "3rd", "4th"), 
                               colnames(datam))
    object@higher <- momentsm
    ## Compute skewness and kurtosis ##
    ## skewness and kurtosis are 1 x r vectors ## 
    skewm <- momentsm[3, ]/momentsm[2, ]^1.5
    kurtm <- momentsm[4, ]/momentsm[2, ]^2
    names(skewm) <- colnames(datam)
    names(kurtm) <- colnames(datam)
    object@skewness <- skewm
    object@kurtosis <- kurtm
    ## Compute corr matrix in case of r > 1 ##
    ## corr is a r x r matrix ##
    if(object@fdata@r  > 1) {
        object@corr <- cor(datam)			
    } else {
        object@corr <- matrix()
    }
    return(object)
}
