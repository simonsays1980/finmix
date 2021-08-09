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

.mcmcoutputhier <- setClass("mcmcoutputhier",
                            representation(hyper = "list"),
                            contains = c("mcmcoutputbase"),
                            validity = function(object) 
                            {
                                ## else: OK
                                TRUE
                            },
                            prototype(hyper = list())
)

setMethod("show", "mcmcoutputhier", 
          function(object)
          {
              cat("Object 'mcmcoutput'\n")
              cat("     class       :", class(object), "\n")
              cat("     M           :", object@M, "\n")
              cat("     burnin      :", object@burnin, "\n")
              cat("     ranperm     :", object@ranperm, "\n")
              cat("     par         : List of", 
                  length(object@par), "\n")
              cat("     log         : List of", 
                  length(object@log), "\n")
              cat("     hyper       : List of",
                  length(object@hyper), "\n")
              cat("     ST          :", 
                  paste(dim(object@ST), collapse = "x"), "\n")
              if (!all(is.na(object@S))) {
                  cat("     S           :", 
                      paste(dim(object@S), collapse = "x"), "\n")
              }
              cat("     NK          :",
                  paste(dim(object@NK), collapse = "x"), "\n")
              cat("     clust       :",
                  paste(dim(object@clust), collapse = "x"), "\n")
              cat("     model       : Object of class", 
                  class(object@model), "\n")
              cat("     prior       : Object of class",
                  class(object@prior), "\n")
          }
)

setMethod( "plotTraces", signature( x     = "mcmcoutputhier", 
                                    dev   = "ANY",
                                    lik   = "ANY",
                                    col   = "ANY" ), 
          function( x, dev = TRUE, lik = 1, col = FALSE, ... ) 
          {
              dist <- x@model@dist
              if ( lik %in% c( 0, 1 ) ) {
                  if ( dist == "poisson" || dist == "cond.poisson" ) {
                      .traces.Poisson.Base.Hier( x, dev )
                  } else if ( dist == "binomial" ) {
                      .traces.Binomial.Base( x, dev )
                  } else if ( dist == "exponential" ) {
                      .traces.Exponential.Base( x, dev )
                  } else if ( dist == "normal" ) {
                      .traces.Normal.Hier( x, dev )
                      .traces.Weights.Base( x, dev, col )
                  } else if ( dist == "student" ) {
                      .traces.Student.Hier( x, dev )
                      .traces.Weights.Base( x, dev, col )
                  } else if ( dist == "normult" ) {
                      .traces.Normult.Hier( x, dev, col )
                      .traces.Weights.Base( x, dev, col )
                  } else if ( dist == "studmult" ) {
                      .traces.Studmult.Hier( x, dev, col )
                      .traces.Weights.Base( x, dev, col )
                  }
              }
              if ( lik %in% c( 1, 2 ) ) {
                  ## log ##
                  .traces.Log.Base( x, dev ) 
              }
          }
)

setMethod("plotHist", signature(x   = "mcmcoutputhier", 
                                dev = "ANY"), 
          function(x, dev = TRUE, ...) 
          {
              dist <- x@model@dist
              if (dist == "poisson") {
                  .hist.Poisson.Base.Hier(x, dev)
              }	else if (dist == "binomial") {
                  .hist.Binomial.Base(x, dev)  
              }
          }
)

setMethod("plotDens", signature(x   = "mcmcoutputhier",
                                dev = "ANY"),
          function(x, dev = TRUE, ...)
          {
              dist <- x@model@dist
              if (dist == "poisson") {
                  .dens.Poisson.Base.Hier(x, dev)
              } else if (dist == "binomial") {
                  .dens.Binomial.Base(x, dev)
              }
          }
)

setMethod("plotPointProc", signature(x      = "mcmcoutputhier",
                                     dev    = "ANY"),
          function(x, dev = TRUE, ...)
          {
              ## Call 'plotPointProc()' from 'mcmcoutputbase'
              callNextMethod(x, dev, ...)
          }
)

setMethod("plotSampRep", signature(x    = "mcmcoutputhier",
                                   dev  = "ANY"),
          function(x, dev = TRUE, ...) 
          {
              ## Call 'plotSampRep()' from 'mcmcoutputbase'
              callNextMethod(x, dev, ...)
          }
)

setMethod("plotPostDens", signature(x   = "mcmcoutputhier",
                                    dev = "ANY"),
          function(x, dev = TRUE, ...)
          {
              ## Call 'plotPostDens()' from 'mcmcoutputbase'
              callNextMethod(x, dev, ...)
          }
)

setMethod( "subseq", signature( object    = "mcmcoutputhier", 
                                index     = "array" ), 
           function( object, index ) 
           {
               ## Call 'subseq()' method from 'mcmcoutputfixhier'
               as( object, "mcmcoutputbase" ) <- callNextMethod( object, index )             
               dist <- object@model@dist
               if ( dist == "poisson" ) {
                   .subseq.Poisson.Hier( object, index )
               } else if ( dist %in% c( "normal", "student" ) ) {
                   .subseq.Norstud.Hier( object, index )
               } else if ( dist %in% c( "normult", "studmult" ) ) {
                   .subseq.Normultstud.Hier( object, index )  
               }
           }
)

setMethod( "swapElements", signature( object  = "mcmcoutputhier", 
                                      index   = "array" ),
          function( object, index ) 
          {
              ## Check arguments, TODO: .validObject ##
              ## Call method 'swapElements()' from 'mcmcoutputbase' 
              callNextMethod( object, index )
          }
)

setMethod( "getHyper", "mcmcoutputhier",
           function( object ) {
               return( object@hyper )
           }
)


## No setters as users are not intended to manipulate ##
## this object ##

### Private functions.
### These functions are not exported.

### Plot
### Plot traces
### Plot traces Poisson: Plots the traces of the MCMC sample
### for the Poisson parameters, the weights and the hyper-
### parameter 'b'.
".traces.Poisson.Base.Hier" <- function(x, dev)
{
    K <- x@model@K
    trace.n <- K * 2
    if (.check.grDevice() && dev) {
        dev.new(title = "Traceplots")
    }
    par(mfrow = c(trace.n, 1), mar = c(1, 0, 0, 0),
        oma = c(4, 5, 4, 4))
    lambda <- x@par$lambda
    for (k in 1:K) {
        plot(lambda[, k], type = "l", axes = F, 
             col = "gray20", xlab = "", ylab = "")
        axis(2, las = 2, cex.axis = 0.7)
        mtext(side = 2, las = 2, bquote(lambda[k = .(k)]),
              cex = 0.6, line = 3)
    }
    weight <- x@weight
    for (k in 1:(K - 1)) {
        plot(weight[, k], type = "l", axes = F, 
             col = "gray47", xlab = "", ylab = "")
        axis(2, las = 2, cex.axis = 0.7)
        mtext(side = 2, las = 2, bquote(eta[k = .(k)]),
              cex = 0.6, line = 3)
    }
    b <- x@hyper$b
    plot(b, type = "l", axes = F,
         col = "gray68", xlab = "", ylab = "")
    axis(2, las = 2, cex.axis = 0.7)
    mtext(side = 2, las = 2, "b", cex = 0.6, line = 3)
    axis(1)
    mtext(side = 1, "Iterations", cex = 0.7, line = 3)	
}

### Plot Histograms
### Plot Histograms Poisson: Plots histograms for
### the Poisson parameters the weights and the hyper-
### parameter b.
".hist.Poisson.Base.Hier" <- function(x, dev)
{
    K <- x@model@K 
    if (.check.grDevice() && dev) {
        dev.new(title = "Histograms")
    }
    lambda <- x@par$lambda
    weight <- x@weight
    b <- x@hyper$b
    vars <- cbind(lambda, weight[, seq(1, K - 1)], b)
    lab.names <- vector("list", 2 * K)
    for (k in seq(1, K)) {
        lab.names[[k]] <- bquote(lambda[.(k)])
    }
    for (k in seq(K + 1, 2 * K - 1)) {
        lab.names[[k]] <- bquote(eta[.(k - K)])
    }
    lab.names[[2 * K]] <- "b"
    .symmetric.Hist(vars, lab.names)
}	

### Plot Densities
### Plot Densities Poisson: Plots Kernel densities for
### the Poisson parameters the weights and the hyper-
### parameter b.
".dens.Poisson.Base.Hier" <- function(x, dev)
{
    K   <- x@model@K
    if (.check.grDevice() && dev) {
        dev.new(title = "Densities")
    }
    lambda  <- x@par$lambda
    weight  <- x@weight
    b       <- x@hyper$b
    vars    <- cbind(lambda, weight[, seq(1, K - 1)], b)
    lab.names <- vector("list", 2 * K)
    for (k in seq(1, K)) {
        lab.names[[k]] <- bquote(lambda[.(k)])
    }
    for (k in seq(K + 1, 2 * K - 1)) {
        lab.names[[k]] <- bquote(eta[.(k - K)])
    }
    lab.names[[2 * K]] <- "b"
    .symmetric.Dens(vars, lab.names)
}

