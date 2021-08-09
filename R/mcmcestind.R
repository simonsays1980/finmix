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

.mcmcestind <- setClass("mcmcestind",
                        representation(eavg = "list"),
                        contains = c("mcmcestfix"),
                        validity = function(object) 
                        {
                            ## else: OK
                            TRUE
                        },
                        prototype(eavg  = list())                                  
)

setClassUnion("mcmcest", 
              c("mcmcestfix",
                "mcmcestind")
)

setMethod("show", "mcmcestind",
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
              cat("     eavg        : List of",
                  length(object@eavg), "\n")
              cat("     sdpost      : List of",
                  length(object@sdpost), "\n")
          }
)

setMethod("Summary", "mcmcestind",
          function(x, ..., na.rm = FALSE)
          {
              dopt  <- getOption("digits")
              options(digits = 4)
              obj   <- x
              K     <- obj@K
              rnames    <- .rownames.Mcmcestind(obj)
              cnames    <- c("Estimates", "Std. Error")
              cat("\n")
              cat("Call: mcmcestimate\n")
              cat("\n")
              if (obj@ranperm) {
                  cat("Method: Random Permutation Gibbs Sampling\n")
              } else {
                  cat("Method: Gibbs Sampling\n")
              }
              cat("\n")
              cat(paste("Number of Iterations: ", obj@M, "\n", sep = ""))
              cat(paste("Number of Burnin Iterations: ", obj@burnin, 
                        "\n", sep = ""))
              cat(paste("Relabeling algorithm used: ", obj@relabel, "\n", 
                        sep = ""))
              cat("\n")
              cat("Parameters:\n")
              cat("\n")
              cat(paste("Component Parameters: ", 
                        .parnames.Mcmcestfix(obj), "\n", sep = ""))
              cat("Weights: eta\n")
              ## MAP ##
              cat("Maximum A Posterior (MAP)\n")
              parout    <- .pars.map.Mcmcestind(obj)
              rownames(parout)  <- rnames
              colnames(parout)  <- cnames              
              print(parout)
              cat("\n")
              cat(paste("Log likelihood: ", sprintf("%.4f", obj@map$log), "\n", sep = ""))
              cat("---\n")
              ## BML ##
              cat("Bayesian Maximum Likelihood (BML)\n")
              parout    <- .pars.bml.Mcmcestind(obj)
              rownames(parout)  <- rnames
              colnames(parout)  <- cnames
              print(parout)
              cat("\n")
              cat(paste("Log likelihood: ", sprintf("%.4f", obj@bml$log), "\n", sep = ""))
              cat("---\n")
              ## IEAVG ##
              cat("Identified Ergodic Average (IEAVG)\n")
              parout    <- .pars.ieavg.Mcmcestind(obj)
              rownames(parout)  <- rnames
              colnames(parout)  <- cnames
              print(parout)
              cat("---\n")
              ## EAVG ##
              cat("Ergodic Average (EAVG)\n")
              parout    <- .pars.eavg.Mcmcestind(obj)
              rownames(parout)  <- rnames
              colnames(parout)  <- cnames
              print(parout)
              cat("---\n")
              options(digits = dopt)
          }
)

## Getters ##
setMethod("getEavg", "mcmcestind", 
          function(object) 
          {
              return(object@eavg)
          }
)

## No setters as users are not intended to manipulate 
## this object.

### Private functions.
### These functions are not exported.

### Summary
### Summary Map estimates: Creates a matrix with Map
### estimates.
".pars.map.Mcmcestind" <- function(obj)
{
    if (obj@dist == "poisson") {
        .pars.map.poisson.Mcmcestind(obj)
    }
}

### Summary Map estimates Poisson: Creates a matrix
### with Map estimates for Poisson parameters.
".pars.map.poisson.Mcmcestind" <- function(obj)
{
    K   <- obj@K
    parout <- matrix(0, nrow = 2 * K, ncol = 2)
    for (k in seq(1, K)) {
        parout[k, 1]    <- obj@map$par$lambda[k]
        parout[k, 2]    <- obj@sdpost$identified$par$lambda[k]
    }
    for (k in seq(1, K)) {
        parout[k + K, 1] <- obj@map$weight[k]
        parout[k + K, 2] <- obj@sdpost$identified$weight[k]
    }
    return(parout)
}

### Summary Bml estimates: Creates a matrix with Bml
### estimates.
".pars.bml.Mcmcestind" <- function(obj)
{
    if (obj@dist == "poisson") {
        .pars.bml.poisson.Mcmcestind(obj)
    }
}

### Summary Bml estimates Poisson: Creates a matrix
### with Bml estimates for Poisson parameters.
".pars.bml.poisson.Mcmcestind" <- function(obj)
{
    K   <- obj@K
    parout <- matrix(0, nrow = 2 * K, ncol = 2)
    for (k in seq(1, K)) {
        parout[k, 1]    <- obj@bml$par$lambda[k]
        parout[k, 2]    <- obj@sdpost$identified$par$lambda[k]
    }
    for (k in seq(1, K)) {
        parout[k + K, 1] <- obj@bml$weight[k]
        parout[k + K, 2] <- obj@sdpost$identified$weight[k]
    }
    return(parout)
}

### Summary Ieavg estimates: Creates a matrix with Ieavg
### estimates.
".pars.ieavg.Mcmcestind" <- function(obj)
{
    if (obj@dist == "poisson") {
        .pars.ieavg.poisson.Mcmcestind(obj)
    }
}

### Summary Bml estimates Poisson: Creates a matrix
### with Bml estimates for Poisson parameters.
".pars.ieavg.poisson.Mcmcestind" <- function(obj)
{
    K   <- obj@K
    parout <- matrix(0, nrow = 2 * K, ncol = 2)
    for (k in seq(1, K)) {
        parout[k, 1]    <- obj@ieavg$par$lambda[k]
        parout[k, 2]    <- obj@sdpost$identified$par$lambda[k]
    }
    for (k in seq(1, K)) {
        parout[k + K, 1] <- obj@ieavg$weight[k]
        parout[k + K, 2] <- obj@sdpost$identified$weight[k]
    }
    return(parout)
}

### Summary Eavg estimates: Creates a matrix with Eavg
### estimates.
".pars.eavg.Mcmcestind" <- function(obj)
{
    if (obj@dist == "poisson") {
        .pars.eavg.poisson.Mcmcestind(obj)
    }
}

### Summary Bml estimates Poisson: Creates a matrix
### with Bml estimates for Poisson parameters.
".pars.eavg.poisson.Mcmcestind" <- function(obj)
{
    K   <- obj@K
    parout <- matrix(0, nrow = 2 * K, ncol = 2)
    for (k in seq(1, K)) {
        parout[k, 1]    <- obj@eavg$par$lambda[k]
        parout[k, 2]    <- obj@sdpost$unidentified$par$lambda[k]
    }
    for (k in seq(1, K)) {
        parout[k + K, 1] <- obj@eavg$weight[k]
        parout[k + K, 2] <- obj@sdpost$unidentified$weight[k]
    }

    return(parout)
}

### Summary rownames: Creates row names for all 
### parameters.
".rownames.Mcmcestind" <- function(obj)
{
    if (obj@dist == "poisson") {
        .rownames.poisson.Mcmcestind(obj)
    }
}

### Summary rownames: Creates row names for 
### each model. 
".rownames.poisson.Mcmcestind" <- function(obj)
{
    rnames  <- rep("", 2 * obj@K)
    for (k in seq(1, obj@K)) {
        rnames[k]  <- paste("lambda ", k, sep = "")
    }
    for(k in seq(obj@K + 1, 2 * obj@K)) {
        rnames[k]  <- paste("eta ", k - obj@K, sep = "")
    }
    return(rnames)
}
