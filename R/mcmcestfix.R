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

.mcmcestfix <- setClass("mcmcestfix",
                        representation(dist        = "character",
                                       K           = "integer",
                                       indicmod    = "character",
                                       burnin      = "integer",
                                       M           = "integer",
                                       ranperm     = "logical",
                                       relabel     = "character",
                                       map         = "list",
                                       bml         = "list",
                                       ieavg       = "list",
                                       sdpost      = "list"),                                       
                        validity = function(object) 
                        {
                            ## else: OK
                            TRUE
                        },
                        prototype(dist      = character(),
                                  K         = integer(),
                                  indicmod  = character(),
                                  burnin    = integer(),
                                  M         = integer(),
                                  ranperm   = logical(),
                                  relabel   = character(),
                                  map       = list(),
                                  bml       = list(),
                                  ieavg     = list(),
                                  sdpost    = list()
                                  )
)

setMethod("show", "mcmcestfix", 
          function(object) 
          {
              cat("Object 'mcmcest'\n")
              cat("     dist        :", object@dist, "\n")
              cat("     K           :", object@K, "\n")
              cat("     indicmod    :", object@indicmod, 
                  "\n")
              cat("     M           :", object@M, "\n")
              cat("     burnin      :", object@burnin, "\n")
              cat("     ranperm     :", object@ranperm, "\n")
              cat("     relabel     :", object@relabel, "\n")
              cat("     map         : List of", 
                  length(object@map), "\n")
              cat("     bml         : List of",
                  length(object@bml), "\n")
              cat("     ieavg       : List of", 
                  length(object@ieavg), "\n")
              cat("     sdpost      : List of",
                  length(object@sdpost), "\n")
          }
)

setMethod("Summary", "mcmcestfix",
          function(x, ..., na.rm = FALSE)
          {
              dopt  <- getOption("digits")
              obj   <- x
              K     <- obj@K
              rnames    <- .rownames.Mcmcestfix(obj)
              cnames    <- c("Estimates", "Std. Error")
              cat("\n")
              cat("Call: mcmcestimate\n")
              cat("\n")
              cat("Method: Gibbs Sampling with fixed indicators\n")
              cat("\n")
              cat(paste("Number of Iterations: ", obj@M, "\n", sep = ""))
              cat(paste("Number of Burnin Iterations: ", obj@burnin, 
                        "\n", sep = ""))
              cat("\n")
              cat("Parameters:\n")
              cat("\n")
              cat(paste("Component Parameters: ", 
                        .parnames.Mcmcestfix(obj), "\n", sep = ""))
              ## MAP ##
              cat("Maximum A Posterior (MAP)\n")
              parout    <- .pars.map.Mcmcestfix(obj)
              rownames(parout)  <- rnames
              colnames(parout)  <- cnames              
              print(parout)
              cat("\n")
              cat(paste("Log likelihood: ", sprintf("%.4f", obj@map$log), "\n", sep = ""))
              cat("---\n")
              ## BML ##
              cat("Bayesian Maximum Likelihood (BML)\n")
              parout    <- .pars.bml.Mcmcestfix(obj)
              rownames(parout)  <- rnames
              colnames(parout)  <- cnames
              print(parout)
              cat("\n")
              cat(paste("Log likelihood: ", sprintf("%.4f", obj@bml$log), "\n", sep = ""))
              cat("---\n")
              ## IEAVG ##
              cat("Identified Ergodic Average (IEAVG)\n")
              parout    <- .pars.ieavg.Mcmcestfix(obj)
              rownames(parout)  <- rnames
              colnames(parout)  <- cnames
              print(parout)
              cat("---\n")
              options(digits = dopt)
          }
)
            
## Getters ##
setMethod("getDist", "mcmcestfix", 
          function(object) 
          {
              return(object@dist)
          }
)

setMethod("getK", "mcmcestfix", 
          function(object) 
          {
              return(object@K)
          }
)

setMethod("getIndicmod", "mcmcestfix", 
          function(object) 
          {
              return(object@indicmod)
          }
)

setMethod("getBurnin", "mcmcestfix",
          function(object)
          {
              return(object@burnin)
          }
)

setMethod("getM", "mcmcestfix",
          function(object) 
          {
              return(object@M)
          }
)

setMethod("getRanperm", "mcmcestfix",
          function(object) 
          {
              return(object)
          }
)

setMethod("getRelabel", "mcmcestfix",
          function(object) 
          {
              return(object@relabel)
          }
)

setMethod("getMap", "mcmcestfix", 
           function(object) 
           {
               return(object@map)
           }
)

setMethod("getBml", "mcmcestfix",
          function(object) 
          {              
              return(object@bml)
          }
)

setMethod("getIeavg", "mcmcestfix", 
          function(object) 
          {
              return(object@ieavg)
          }
)

setMethod("getSdpost", "mcmcestfix", 
          function(object) 
          {
              return(object@sdpost)
          }
)

## No setters as users are not intended to manipulate
## this object

### Private functions.
### These functions are not exported.

### Summary
### Summary Map estimates: Creates a matrix with Map
### estimates.
".pars.map.Mcmcestfix" <- function(obj)
{
    if (obj@dist == "poisson") {
        .pars.map.poisson.Mcmcestfix(obj)
    }
}

### Summary Map estimates Poisson: Creates a matrix
### with Map estimates for Poisson parameters.
".pars.map.poisson.Mcmcestfix" <- function(obj)
{
    parout <- matrix(0, nrow = obj@K, ncol = 2)
    for (k in seq(1, obj@K)) {
        parout[k, 1]    <- obj@map$par$lambda[k]
        parout[k, 2]    <- obj@sdpost$identified$par$lambda[k]
    }
    return(parout)
}

### Summary Bml estimates: Creates a matrix with Bml
### estimates.
".pars.bml.Mcmcestfix" <- function(obj)
{
    if (obj@dist == "poisson") {
        .pars.bml.poisson.Mcmcestfix(obj)
    }
}

### Summary Bml estimates Poisson: Creates a matrix
### with Bml estimates for Poisson parameters.
".pars.bml.poisson.Mcmcestfix" <- function(obj)
{
    parout <- matrix(0, nrow = obj@K, ncol = 2)
    for (k in seq(1, obj@K)) {
        parout[k, 1]    <- obj@bml$par$lambda[k]
        parout[k, 2]    <- obj@sdpost$identified$par$lambda[k]
    }
    return(parout)
}

### Summary Ieavg estimates: Creates a matrix with Ieavg
### estimates.
".pars.ieavg.Mcmcestfix" <- function(obj)
{
    if (obj@dist == "poisson") {
        .pars.ieavg.poisson.Mcmcestfix(obj)
    }
}

### Summary Bml estimates Poisson: Creates a matrix
### with Bml estimates for Poisson parameters.
".pars.ieavg.poisson.Mcmcestfix" <- function(obj)
{
    parout <- matrix(0, nrow = obj@K, ncol = 2)
    for (k in seq(1, obj@K)) {
        parout[k, 1]    <- obj@ieavg$par$lambda[k]
        parout[k, 2]    <- obj@sdpost$identified$par$lambda[k]
    }
    return(parout)
}

### Summary rownames: Creates rownames for the summary.
".rownames.Mcmcestfix" <- function(obj)
{
    if (obj@dist == "poisson") {
        .rownames.poisson.Mcmcestfix(obj)
    }
}

### Summary rownames Poisson: Creates the row names
### for the summary of Poisson estimates.
".rownames.poisson.Mcmcestfix" <- function(obj)
{    
    rnames <- rep("", obj@K)
    for (k in seq(1, obj@K)) {
        rnames[k] <- paste("lambda ", k, sep = "")
    }
    return(rnames)
}

### Summary parameter names: Creates parameter
### names for the components.
".parnames.Mcmcestfix" <- function(obj) 
{
    if (obj@dist == "poisson") {
        parnames <- c("lambda")
    } 
    return(parnames)
}
