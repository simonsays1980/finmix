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

.studmultmodelmoments <- setClass("studmultmodelmoments", 
                                  representation(B      = "array",
                                                 W      = "array",
                                                 Rdet   = "numeric",
                                                 Rtr    = "numeric",
                                                 corr   = "array"
                                                     ),
                                  contains = c("cmodelmoments"),
                                  validity = function(object) {
                                      ## else: OK 
                                      TRUE
                                  },
                                  prototype(
                                            B   = array(),
                                            W   = array(),
                                            Rdet= numeric(),
                                            Rtr = numeric(),
                                            corr= array()
                                            )
)

setMethod("initialize", "studmultmodelmoments",
          function(.Object, ..., model) 
          {
              .Object <- callNextMethod(.Object, ..., model = model)
              generateMoments(.Object)
          }
)

setMethod("generateMoments", "studmultmodelmoments",
          function(object) 
          {
              .generateMomentsStudmult(object)
          }
)

setMethod("show", "studmultmodelmoments",
          function(object) 
          {
              cat("Object 'modelmoments'\n")
              cat("     mean        : Vector of", 
                  length(object@mean), "\n")
              cat("     var         :",
                  paste(dim(object@var), collapse = "x"), "\n")
              cat("     higher      :", 
                  paste(dim(object@higher), collapse = "x"),"\n")
              cat("     skewness    : Vector of",
                  length(object@skewness), "\n")
              cat("     kurtosis    : Vector of",
                  length(object@kurtosis), "\n")
              cat("     B           :", 
                  paste(dim(object@B), collapse = "x"), "\n")
              cat("     W           :",
                  paste(dim(object@W), collapse = "x"), "\n")
              cat("     Rdet        :", object@Rdet, "\n")
              cat("     Rtr         :", object@Rtr, "\n")
              cat("     corr        :", 
                  paste(dim(object@corr), collapse = "x"), "\n")
              cat("     model       : Object of class",
                  class(object@model), "\n")
          }
)

## Getters ##
setMethod("getB", "studmultmodelmoments",
          function(object) 
          {
              return(object@B)
          }
)

setMethod("getW", "studmultmodelmoments",
          function(object) 
          {
              return(object@W)
          }
)

setMethod("getRdet", "studmultmodelmoments",
          function(object) 
          {
              return(object@Rdet)
          }
)

setMethod("getRtr", "studmultmodelmoments",
          function(object) 
          {
              return(object@Rtr)
          }
)

setMethod("getCorr", "studmultmodelmoments",
          function(object) 
          {
              return(object@corr)
          }
)

## No setters as users are not intended to manipulate ##
## this object ##

### Private functions
### These function are not exported
".generateMomentsStudmult" <- function(object) 
{
    mu          <- object@model@par$mu
    sigma       <- object@model@par$sigma
    df          <- object@model@par$df
    weight      <- object@model@weight
    names       <- rep("", object@model@r)
    for (i in seq(1, object@model@r)) {
        names[i] <- paste("r=", i, sep = "")
    }
    object@mean <- apply(apply(mu, 1, '*', weight)
                         , 2, sum, na.rm = TRUE)
    if(all(df > 2)) {
        object@W    <- apply(sweep(sigma, MARGIN = 3, 
                                   weight * df/(df - 2), '*'),
                             c(1, 2), sum, na.rm = TRUE)
        object@var  <- object@W + apply(apply(mu, 2, 
                                              tcrossprod, mu),
                                        1, '*',
                                        weight)
        object@var  <- object@var - object@mean %*% t(object@mean)
        diffm       <- mu - object@mean
        object@B    <- apply(apply(diffm, 1, tcrossprod, diffm),
                             1, '*', weight)
        cd          <- diag(1/diag(object@var)^.5)
        object@corr <- cd %*% object@var %*% cd
        object@Rtr  <- 1 - sum(diag(object@W))/sum(diag(object@var))
        object@Rdet <- 1 - det(object@W)/det(object@var)
    } else {
        r <- object@model@r
        object@W        <- array(NaN, dim = c(r, r)) 
        object@var      <- array(NaN, dim = c(r, r))
        object@B        <- array(NaN, dim = c(r, r))
        object@Rdet     <- NaN
        object@Rtr      <- NaN
        object@corr     <- array(NaN, dim = c(r, r))
    }
    names(object@mean)      <- names
    colnames(object@var)    <- names
    rownames(object@var)    <- names
    colnames(object@B)      <- names
    rownames(object@B)      <- names
    colnames(object@W)      <- names
    rownames(object@W)      <- names
    colnames(object@corr)   <- names
    rownames(object@corr)   <- names
    highm           <- array(0, dim = c(4, object@model@r))
    dimnames(highm) <- list(c("1st", "2nd", "3rd", "4th"), names)
    for (i in seq(1, object@model@r)) {
        marmodel    <- mixturemar(object@model, i)
        highm[, i]  <- .mixturemoments.student(marmodel, 4, 
                                               object@mean[i]) 
    }
    object@higher   <- highm
    object@skewness <- object@higher[3, ]/object@higher[2, ]^1.5
    object@kurtosis <- object@higher[4, ]/object@higher[2, ]^2
    return(object)
}


    
