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

.csdatamoments <- setClass("csdatamoments",
                           representation(B       = "vector",
                                          W       = "vector",
                                          T       = "vector",
                                          R       = "numeric",
                                          Rtr     = "numeric",
                                          Rdet    = "numeric"),
                           contains = c("sdatamoments"),
                           validity = function(object) 
                           {
                               ## else: ok
                               TRUE
                           },
                           prototype(B    = vector("numeric"),
                                     W    = vector("numeric"),
                                     T    = vector("numeric"),
                                     R    = numeric(),
                                     Rtr  = numeric(),
                                     Rdet = numeric()
                                     )
)

setClassUnion("csdatamomentsOrNULL", members = c("csdatamoments", "NULL"))

setMethod("initialize", "csdatamoments",
          function(.Object, ..., value = fdata()) 
          {
              .Object <- callNextMethod(.Object, ..., value = value)
              if(hasY(value) && hasS(value)) {
                  .Object <- generateMoments(.Object)
              }
              return(.Object)
          }
)

setMethod("generateMoments", "csdatamoments",
          function(object) 
          {
              .generateCsdatamoments(object)
          }
)

setMethod("show", "csdatamoments", 
          function(object) 
          {
              cat("Object 'sdatamoments'\n")
              cat("     B           : Vector of", 
                  length(object@B), "\n")
              cat("     W           : Vector of",
                  length(object@W), "\n")
              cat("     T           : Vector of",
                  length(object@T), "\n")
              if (object@fdata@r > 1) {
                  cat("     Rdet        :", object@Rdet, "\n")
                  cat("     Rtr         :", object@Rtr, "\n")
              }
              cat("     gmoments    : Object of class", 
                  class(object@gmoments), "\n")
              cat("     fdata        : Object of class", 
                  class(object@fdata), "\n")
          }
)

## Getters ##
setMethod("getGmoments", "csdatamoments", 
          function(object) 
          {
              return(object@gmoments)
          }
)

setMethod("getWK", "csdatamoments", 
          function(object) 
          {
              return(object@WK)
          }
)

setMethod("getVar", "csdatamoments", 
          function(object) 
          {
              return(object@var)
          }
)

setMethod("getB", "csdatamoments", 
          function(object) 
          {
              return(object@B)
          }
)

setMethod("getW", "csdatamoments", 
          function(object) 
          {
              return(object@W)
          }
)

setMethod("getT", "csdatamoments", 
          function(object) 
          {
              return(object@T)
          }
)

setMethod("getR", "csdatamoments",
          function(object) 
          {
              return(object@R)
          }
)

setMethod("getRtr", "csdatamoments", 
          function(object) 
          {
              return(object@Rtr)
          }
)

setMethod("getRdet", "csdatamoments", 
          function(object) 
          {
              return(object@Rdet)
          }
)

setMethod("getFdata", "csdatamoments", 
          function(object) 
          {
              return(object@fdata)				
          }
)

## Setters ##
## No setters, as it users are adviced not to manipulate moment objects ##

### Private functions
### These functions are not exported
".generateCsdatamoments" <- function(object)
{
    ## enforce column.wise ordering ##
    datam   <- getColY(object@fdata)
    classm  <- getColS(object@fdata)
    ## Calculate the between-group variance ##
    ## 'B' is an r x r matrix ##
    gmeans <- object@gmoments@mean
    nkm <- object@gmoments@NK
    ## Calculate the total heterogeneity ##
    ## 'T' is an r x r array ##
    object@T <- var(datam, na.rm = TRUE) * nrow(datam)
    ## Calculate the within-group heterogeneity ##
    ## 'W' is an r x r array ##
    wkm <- object@gmoments@WK
    object@W <- apply(wkm, c(1, 2), sum, na.rm = TRUE)
    ## Calculate between-group heterogeneity ##
    ## 'B' is an r x r array ##
    object@B <- object@T - object@W
    ## Calculate coefficient of determination ##
    ## 'Rtr' is an 1 x 1 numeric ##
    ## 'Rdet' is an 1 x 1 numeric ##
    if (object@data@r > 1) {
        r <- NA
        object@R <- as.numeric(r)
        object@Rtr <- 1 - sum(diag(object@W), na.rm = TRUE) / 
        sum(diag(object@T), na.rm = TRUE) 
        object@Rdet <- 1 - det(object@W)/det(object@T)
    } else {
        rtr <- NA
        rdet <- NA
        object@Rtr <- as.numeric(rtr)
        object@Rdet <- as.numeric(rdet)
        object@R <- 1 - object@W[1]/object@T[1] 
    }
    return(object)
}
