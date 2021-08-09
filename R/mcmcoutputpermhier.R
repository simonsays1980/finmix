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

.mcmcoutputpermhier <- setClass("mcmcoutputpermhier",
                                contains = c("mcmcpermind", 
                                             "mcmcoutputhier"),
                                validity = function(object) 
                                {
                                    ## else: OK
                                    TRUE
                                }
)

setMethod("initialize", "mcmcoutputpermhier",
          function(.Object, mcmcoutput, Mperm = integer(), 
                   parperm = list(), relabel = character(), 
                   weightperm = array(), logperm = list(), 
                   entropyperm = array(), STperm = array(), 
                   Sperm = array(), NKperm = array()) 
          {
              .Object@M             <- mcmcoutput@M
              .Object@burnin        <- mcmcoutput@burnin
              .Object@ranperm       <- mcmcoutput@ranperm
              .Object@par           <- mcmcoutput@par
              .Object@weight        <- mcmcoutput@weight
              .Object@log           <- mcmcoutput@log
              .Object@hyper         <- mcmcoutput@hyper
              .Object@ST            <- mcmcoutput@ST
              .Object@S             <- mcmcoutput@S
              .Object@NK            <- mcmcoutput@NK
              .Object@clust         <- mcmcoutput@clust
              .Object@model         <- mcmcoutput@model
              .Object@prior         <- mcmcoutput@prior
              .Object@Mperm         <- Mperm
              .Object@parperm       <- parperm
              .Object@relabel       <- relabel
              .Object@weightperm    <- weightperm
              .Object@logperm       <- logperm
              .Object@entropyperm   <- entropyperm
              .Object@STperm        <- STperm
              .Object@Sperm         <- Sperm
              .Object@NKperm        <- NKperm
              .Object
          }
)

setMethod("show", "mcmcoutputpermhier", 
          function(object)
          {
              cat("Object 'mcmcoutputperm'\n")
              cat("     class       :", class(object), "\n")
              cat("     M           :", object@M, "\n")
              cat("     burnin      :", object@burnin, "\n")
              cat("     ranperm     :", object@ranperm, "\n")
              cat("     relabel     :", object@relabel, "\n")
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
              cat("     Mperm       :", object@Mperm, "\n")
              cat("     parperm     : List of", 
                  length(object@parperm), "\n")
              cat("     weightperm  :",
                  paste(dim(object@weightperm), collapse = "x"), "\n")
              cat("     logperm     : List of",
                  length(object@logperm), "\n")
              cat("     entropyperm :",
                  paste(dim(object@entropyperm), collapse = "x"), "\n")
              cat("     STperm      :",
                  paste(dim(object@STperm), collapse = "x"), "\n")
              if (!all(is.na(object@Sperm))) {
                  cat("     Sperm       :",
                      paste(dim(object@Sperm), collapse = "x"), "\n")
              }
              cat("     NKperm      :", 
                  paste(dim(object@NKperm), collapse = "x"), "\n")
              cat("     model       : Object of class", 
                  class(object@model), "\n")
              cat("     prior       : Object of class",
                  class(object@prior), "\n")
          }
)

setMethod( "plotTraces", signature( x     = "mcmcoutputpermhier",
                                    dev   = "ANY",
                                    lik   = "ANY",
                                    col   = "ANY" ), 
          function( x, dev = TRUE, lik = 1, col = FALSE, ... ) 
          {
              dist <- x@model@dist
              if ( lik %in% c( 0, 1 ) ) {
                  if ( dist == "poisson" ) {
                      .permtraces.Poisson.Base.Hier( x, dev )
                  } else if (dist == "binomial") {
                      .permtraces.Binomial.Base( x, dev )
                  } else if ( dist == "exponential" ) {
                      .permtraces.Exponential.Base( x, dev )
                  } else if ( dist == "normal" ) {
                      .permtraces.Normal.Hier( x, dev )
                      .permtraces.Weights.Base( x, dev, col )
                  } else if ( dist == "student" ) {
                      .permtraces.Student.Hier( x, dev )
                      .permtraces.Weights.Base( x, dev, col )
                  } else if ( dist == "normult" ) {
                      .permtraces.Normult.Hier( x, dev, col )
                      .permtraces.Weights.Base( x, dev, col )
                  } else if ( dist == "studmult" ) {
                      .permtraces.Studmult.Hier( x, dev, col )
                      .permtraces.Weights.Base( x, dev, col )
                  }
              }	
              if ( lik %in% c( 1, 2 ) ) {
                  ## log ##
                  .permtraces.Log.Base( x, dev )
              }
          }
)

setMethod("plotHist", signature(x = "mcmcoutputpermhier", 
                                dev = "ANY"), 
          function(x, dev = TRUE, ...) 
          {
              dist <- x@model@dist
              if(dist == "poisson") {
                  .permhist.Poisson.Base.Hier(x, dev)
              }	else if (dist == "binomial") {
                  .permhist.Binomial.Base(x, dev)
              }
          }
)

setMethod("plotDens", signature(x = "mcmcoutputpermhier", 
                                dev = "ANY"), 
          function(x, dev = TRUE, ...) 
          {
              dist <- x@model@dist
              if (dist == "poisson") {
                  .permdens.Poisson.Base.Hier(x, dev)
              }	else if (dist == "binomial") {
                  .permhist.Binomial.Base(x, dev)
              }
           }
)

setMethod("plotPointProc", signature(x      = "mcmcoutputpermhier",
                                     dev    = "ANY"),
          function(x, dev = TRUE, ...)
          {
              dist <- x@model@dist
              if (dist == "poisson") {
                  .permpointproc.Poisson(x, dev)
              } else if (dist == "binomial") {
                  .permpointproc.Binomial(x, dev)
              }
          }
)

setMethod("plotSampRep", signature(x    = "mcmcoutputpermhier",
                                   dev  = "ANY"),
          function(x, dev, ...) 
          {
              dist <- x@model@dist
              if (dist == "poisson") {
                  .permsamprep.Poisson(x, dev)
              } else if (dist == "binomial") {
                  .permsamprep.Binomial(x, dev)
              }
          }
)

setMethod("plotPostDens", signature(x   = "mcmcoutputpermhier",
                                    dev = "ANY"),
          function(x, dev = TRUE, ...) 
          {
              dist <- x@model@dist
              if (dist == "poisson") {
                  .permpostdens.Poisson(x, dev)
              } else if (dist == "binomial") {
                  .permpostdens.Binomial(x, dev)
              }
          }
)

### Private functions.
### These functions are not exported.

### Plot
### Traces
### Traces Poisson: Plots the traces for all Poisson 
### parameters, the weights and the hpyer-parameter 'b'.
".permtraces.Poisson.Base.Hier" <- function(x, dev)
{
    K <- x@model@K
		trace.n <- K * 2
    if (.check.grDevice() && dev) {
        dev.new(title = "Traceplots (permuted)")
    }
    par(mfrow = c(trace.n, 1), mar = c(1, 0, 0, 0),
        oma = c(4, 5, 4, 4))
    lambda <- x@parperm$lambda
    for (k in 1:K) {
        plot(lambda[, k], type = "l", axes = F, 
             col = "gray20", xlab = "", ylab = "")
        axis(2, las = 2, cex.axis = 0.7)
        mtext(side = 2, las = 2, bquote(lambda[k = .(k)]),
              cex = 0.6, line = 3)
    }
    weight <- x@weightperm
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

### Histograms
### Histograms Poisson: plots histograms for all Poisson parameters,
### the weights and the hyper-parameter 'b'.
".permhist.Poisson.Base.Hier" <- function(x, dev)
{
    K <- x@model@K 
    if (.check.grDevice() && dev) {
        dev.new(title = "Histograms (permuted)")
    }
    lambda <- x@parperm$lambda
    weight <- x@weightperm
    b <- x@hyper$b
    vars        <- cbind(lambda, weight[, seq(1:(K - 1))], b)
    lab.names   <- vector("list", 2 * K)
    for (k in 1:K) {
        lab.names[[k]] <- bquote(lambda[.(k)])
    }
    for (k in (K + 1):(2 * K - 1)) {
        lab.names[[k]] <- bquote(eta[.(k - K)])
    }
    lab.names[[2 * K]] <- "b"
    .symmetric.Hist(vars, lab.names)
}

### Densities
### Densities Poisson: plots Kernel densities for all Poisson 
### parameters, the weights and the hyper-parameter 'b'.
".permdens.Poisson.Base.Hier" <- function(x, dev)
{
    K <- x@model@K 
    if (.check.grDevice() && dev) {
        dev.new(title = "Histograms (permuted)")
    }
    lambda <- x@parperm$lambda
    weight <- x@weightperm
    b <- x@hyper$b
    vars        <- cbind(lambda, weight[, seq(1:(K - 1))], b)
    lab.names   <- vector("list", 2 * K)
    for (k in 1:K) {
        lab.names[[k]] <- bquote(lambda[.(k)])
    }
    for (k in (K + 1):(2 * K - 1)) {
        lab.names[[k]] <- bquote(eta[.(k - K)])
    }
    lab.names[[2 * K]] <- "b"
    .symmetric.Dens(vars, lab.names)
}
