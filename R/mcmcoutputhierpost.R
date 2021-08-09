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

.mcmcoutputhierpost <- setClass("mcmcoutputhierpost", 
                                representation(post = "list"),
                                contains = c("mcmcoutputhier"),
                                validity = function(object) 
                                {
                                    ## else: OK
                                    TRUE
                                }, 
                                prototype(post  = list())
)

## Set 'mcmcoutput' to the virtual class inheriting ##
## to each other 'mcmcoutput' class. 			    ##
## This is done to simplify dispatching methods.	##
setClassUnion("mcmcoutput", 
	c(
		"mcmcoutputfix",
		"mcmcoutputfixhier",
		"mcmcoutputfixpost",
		"mcmcoutputfixhierpost",
		"mcmcoutputbase",
		"mcmcoutputhier",
		"mcmcoutputpost",
		"mcmcoutputhierpost")
)

setMethod("show", "mcmcoutputhierpost", 
          function(object)
          {
              cat("Object 'mcmcoutput'\n")
              cat("     class       :", class(object), "\n")
              cat("     M           :", object@M, "\n")
              cat("     burnin      :", object@burnin, "\n")
              cat("     ranperm     :", object@ranperm, "\n")
              cat("     par         : List of", 
                  length(object@par), "\n")
              cat("     weight      :",
                  paste(dim(object@weight), collapse = "x"), "\n")
              cat("     log         : List of", 
                  length(object@log), "\n")
              cat("     hyper       : List of",
                  length(object@hyper), "\n")
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
              cat("     model       : Object of class", 
                  class(object@model), "\n")
              cat("     prior       : Object of class",
                  class(object@prior), "\n")
          }
)

setMethod( "plotTraces", signature( x     = "mcmcoutputhierpost", 
                                    dev   = "ANY",
                                    lik   = "ANY",
                                    col   = "ANY" ), 
          function(x, dev = TRUE, lik = 1, col = FALSE, ...) 
          {
              ## Call method 'plot()' from 'mcmcoutputhier'
              callNextMethod(x, dev, lik, col, ...)
          }
)

setMethod("plotHist", signature(x   = "mcmcoutputhierpost", 
                                dev = "ANY"), 
          function(x, dev = TRUE, ...) 
          {
              ## Call method 'plotHist()' from 'mcmcoutputhier'
              callNextMethod(x, dev, ...)
          }
)

setMethod("plotDens", signature(x   = "mcmcoutputhierpost", 
                                dev = "ANY"), 
          function(x, dev = TRUE, ...) 
          {
              ## Call 'plotHist()' from 'mcmcoutputhier'
              callNextMethod(x, dev, ...)
          }
)

setMethod("plotPointProc", signature(x      = "mcmcoutputhierpost",
                                     dev    = "ANY"),
          function(x, dev = TRUE, ...)
          {
              ## Call 'plotPointProc()' from 'mcmcoutputhier'
              callNextMethod(x, dev, ...)
          }
)

setMethod("plotSampRep", signature(x    = "mcmcoutputhierpost",
                                   dev  = "ANY"),
          function(x, dev = TRUE, ...) 
          {
              ## Call 'plotSampRep()' from 'mcmcoutputhier'
              callNextMethod(x, dev, ...)
          }
)

setMethod("plotPostDens", signature(x   = "mcmcoutputhierpost",
                                    dev = "ANY"),
          function(x, dev = TRUE, ...)
          {
              ## Call 'plotPostDens()' from 'mcmcoutputhier'
              callNextMethod(x, dev, ...)
          }
)

setMethod( "subseq", signature( object = "mcmcoutputhierpost", 
                                index  = "array"), 
          function( object, index ) 
          {
              ## Call 'subseq()' method from 'mcmcoutputhier'
              as( object, "mcmcoutputhier" ) <- callNextMethod( object, index )
              # Change owned slots #
              dist      <- object@model@dist
              if ( dist == "poisson" ) {
                  .subseq.Poisson.Post( object, index )
              } else if ( dist == "binomial" ) {
                  .subseq.Binomial.Mcmcoutputfixpost( object, index )
              } else if ( dist %in% c( "normal", "student" ) ) {
                  .subseq.Norstud.Mcmcoutputfixpost( object, index )
              } else if ( dist %in% c( "normult", "studmult" ) ) {
                  .subseq.Normultstud.Mcmcoutputfixpost( object, index ) 
              }

          }
)

setMethod("swapElements", signature( object = "mcmcoutputhierpost", 
                                     index = "array" ),
          function( object, index ) 
          {
              ## Check arguments, TODO: .validObject ##
              if ( object@model@K == 1 ) {
                  return( object )
              } else {
                  dist <- object@model@dist
                  ## Call method 'swapElements()' from 'mcmcoutputhier' 
                  object <- callNextMethod( object, index )
                  if ( dist == "poisson" ) {
                      .swapElements.Poisson.Post( object, index )
                  } else if ( dist == "binomial" ) {
                      .swapElements.binomial.Mcmcoutputfixpost( object, index )
                  } else if ( dist %in% c(  "normal", "student" ) ) {
                      .swapElements.Norstud.Mcmcoutputfixpost( object, index )
                  } else if ( dist %in% c( "normult", "studmult" ) ) {
                      .swapElements.Normultstud( object, index )
                  }
              }
          }
)

setMethod( "getPost", "mcmcoutputhierpost",  
           function( object ) 
           {
               return( object@post )
           }
)

## No setters as users are not intended to manipuate ##
## this object ##
