# Copyright (C) 2013 Lars Simon Zehnder
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

.groupmoments <- setClass("groupmoments", 
                          representation(NK          = "array",
                                         mean        = "matrix",
                                         WK          = "array",
                                         var         = "array",
                                         fdata       = "fdata"), 
                          validity = function(object) {
                              ## else: ok
                              TRUE
                          },
                          prototype(NK       = array(),
                                    mean     = matrix(),
                                    WK       = array(),
                                    var      = array(),
                                    fdata    = fdata()
                                    )
)

"groupmoments" <- function(value = fdata()) 
{
    hasY(value, verbose = TRUE)
    hasS(value, verbose = TRUE)
    .groupmoments(value = value)
}

## initializes by immediately calling method ##
## 'generateMoments' ##
setMethod("initialize", "groupmoments",
          function(.Object, ..., value) 
          {
              .Object@fdata <- value
              generateMoments(.Object)
          }
)

setMethod("generateMoments", "groupmoments",
          function(object) 
          {
              .generateGroupMoments(object)
          }
)

## R usual 'show' function ##
setMethod("show", "groupmoments", 
          function(object) 
          {
              cat("Object 'groupmoments'\n")
              cat("     NK          : Vector of",
                  length(object@NK), "\n")
              cat("     mean        :",
                  paste(dim(object@mean), collapse = "x"), "\n")
              cat("     WK          :",
                  paste(dim(object@WK), collapse = "x"), "\n")
              cat("     var         :",
                  paste(dim(object@var), collapse = "x"), "\n")
              cat("     fdata       : Object of class",
                  class(object@fdata), "\n")
          }
)
 
## R usual Getters ##
setMethod("getNK", "groupmoments", 
          function(object) 
          {
              return(object@NK)					
          }
)

setMethod("getMean", "groupmoments", 
          function(object) 
          {
              return(object@mean)	
          }
)

setMethod("getWK", "groupmoments", 
          function(object) 
          {
              return(object@WK)
          }
)

setMethod("getVar", "groupmoments", 
          function(object) 
          {
              return(object@var)
          }
)

setMethod("getFdata", "groupmoments", 
          function(object) 
          {
              return(object@fdata)		
          }
)
## No setters as user are not intended to manipulate this  ##
## object ##

### Private functions
### These functions are not exported
".generateGroupMoments" <- function(object) 
{
    if(!hasS(object@fdata)) {
        return(object)
    }

    ## Compute group sizes ##
    ## enforce column-wise ordering ##

    datam   <- getColY(object@fdata)
    classm  <- getColS(object@fdata)
    ## Calculate group sizes and group means ##
    ## 'NK' is an 1 x K vector ##
    ## 'groupmean' is an r x K matrix ##
    level.set <- as.numeric(levels(factor(classm)))
    K <- length(level.set)
    r <- ncol(datam)
    comp <- matrix(rep(classm, K), ncol = K) == matrix(seq(1,K), 
                                                       nrow = nrow(datam),
                                                       ncol = K,
                                                       byrow = TRUE)
    names <- rep("", K)
    for (k in seq(1, K)) {
        names[k] <- paste("k=", k, sep = "")
    }
    object@NK <- as.array(apply(comp, 2, sum))
    dimnames(object@NK) <- list(names)
    gmeans <- matrix(NA, nrow = r, ncol = K)
    for (i in seq(1,r)) {
        gmeans[i, ] <- (t(datam[,i]) %*% comp)/t(object@NK)
    }
    colnames(gmeans) <- names
    rownames(gmeans) <- colnames(datam)
    object@mean <- gmeans
    wkm <- array(NA, dim = c(r, r, K))
    varm <- array(NA, dim = c(r, r, K))
    for (k in seq(1, K)) {
        group.demeaned <- (datam - rep(gmeans[,k], each = nrow(datam))) * comp[, k]
        wkm[,, k] <- t(group.demeaned) %*% group.demeaned
        varm[,, k] <- wkm[,, k]/object@NK[k]
    }
    dimnames(wkm) <- list(colnames(datam), colnames(datam), names)
    dimnames(varm) <- list(colnames(datam), colnames(datam), names)
    object@WK <- wkm
    object@var <- varm
    return(object)
}

