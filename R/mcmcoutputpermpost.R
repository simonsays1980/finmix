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

.mcmcoutputpermpost <- setClass("mcmcoutputpermpost",
                                contains = c("mcmcpermindpost", 
                                             "mcmcoutputpost"),
                                validity = function(object) 
                                {
                                    ## else: OK
                                    TRUE
                                }
)

setMethod("initialize", "mcmcoutputpermpost",
          function(.Object, mcmcoutput, Mperm = integer(), 
                   parperm = list(), relabel = character(),
                   weightperm = array(), logperm = list(), 
                   postperm = list(), entropyperm = array(), 
                   STperm = array(), Sperm = array(), 
                   NKperm = array()) 
          {
              .Object@M             <- mcmcoutput@M
              .Object@burnin        <- mcmcoutput@burnin
              .Object@ranperm       <- mcmcoutput@ranperm
              .Object@par           <- mcmcoutput@par
              .Object@weight        <- mcmcoutput@weight
              .Object@log           <- mcmcoutput@log
              .Object@post          <- mcmcoutput@post
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
              .Object@postperm      <- postperm
              .Object@entropyperm   <- entropyperm
              .Object@STperm        <- STperm
              .Object@Sperm         <- Sperm
              .Object@NKperm        <- NKperm
              .Object
          }
)

setMethod("show", "mcmcoutputpermpost", 
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
              cat("     weight      :",
                  paste(dim(object@weight), collapse = "x"), "\n")
              cat("     log         : List of", 
                  length(object@log), "\n")
              cat("     post        : List of",
                  length(object@post), "\n")
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
              cat("     postperm    : List of",
                  length(object@postperm), "\n")
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

setMethod( "plotTraces", signature( x     = "mcmcoutputpermpost", 
                                    dev   = "ANY",
                                    lik   = "ANY",
                                    col   = "ANY" ), 
          function( x, dev = TRUE, lik = 1, col = FALSE, ... ) 
          {
              dist <- x@model@dist
              if ( lik %in% c( 0, 1 ) ) {
                  if ( dist == "poisson" ) {
                      .permtraces.Poisson.Base( x, dev )
                  } else if ( dist == "binomial" ) {
                      .permtraces.Binomial.Base( x, dev )
                  } else if ( dist == "exponential" ) {
                      .permtraces.Exponential.Base( x, dev )
                  } else if ( dist == "normal" ) {
                      .permtraces.Normal( x, dev )
                      .permtraces.Weights.Base( x, dev, col )
                  } else if ( dist == "student" ) {
                      .permtraces.Student( x, dev ) 
                      .permtraces.Weights.Base( x, dev, col )
                  } else if ( dist == "normult" ) {
                      .permtraces.Normult( x, dev, col )
                      .permtraces.Weights.Base( x, dev, col )
                  } else if ( dist == "studmult" ) {
                      .permtraces.Studmult( x, dev, col )
                      .permtraces.Weights.Base(x, dev, col )
                  }

              }
              if ( lik %in% c( 1, 2 ) ) {
                  ## log ##
                  .permtraces.Log.Base( x, dev, col ) 
              }
          }
)

setMethod("plotHist", signature(x   = "mcmcoutputpermpost", 
                                dev = "ANY"), 
          function(x, dev = TRUE, ...) 
          {
              dist <- x@model@dist
              if (dist == "poisson") {
                  .permhist.Poisson.Base(x, dev)
              }	else if (dist == "binomial") {
                  .permhist.Binomial.Base(x, dev)
              }
          }
)

setMethod("plotDens", signature(x   = "mcmcoutputpermpost", 
                                dev = "ANY"), 
          function(x, dev = TRUE, ...) 
          {
              dist <- x@model@dist
              if (dist == "poisson") {
                  .permdens.Poisson.Base(x, dev)
              }	else if (dist == "binomial") {   
                  .permdens.Binomial.Base(x, dev)
              }
          }
)

setMethod("plotPointProc", signature(x      = "mcmcoutputpermpost",
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

setMethod("plotSampRep", signature(x    = "mcmcoutputpermpost",
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

setMethod("plotPostDens", signature(x   = "mcmcoutputpermpost",
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

