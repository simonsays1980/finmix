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

.mcmcoutputfixhierpost <- setClass( "mcmcoutputfixhierpost", 
                                    representation( post = "list" ),
                                    contains = c("mcmcoutputfixhier" ),
                                    validity = function( object ) 
                                    {
                                        ## else: OK
                                        TRUE
                                    },
                                    prototype( post  = list() )
)

setMethod("show", "mcmcoutputfixhierpost", 
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
              cat("     post        : List of",
                  length(object@post), "\n")
              cat("     model       : Object of class",
                  class(object@model), "\n")
              cat("     prior       : Object of class",
                  class(object@prior), "\n")
          }
)

setMethod( "plotTraces", signature( x     = "mcmcoutputfixhierpost", 
                                    dev   = "ANY",
                                    lik   = "ANY",
                                    col   = "ANY" ), 
          function(x, dev = TRUE, lik = 1, col = FALSE, ...) 
          {
              ## Call method 'plot()' from 'mcmcoutputfixhier'
              callNextMethod(x, dev, lik, col, ...)
          }
)

setMethod( "plotHist", signature( x   = "mcmcoutputfixhierpost", 
                                  dev = "ANY" ), 
           function( x, dev = TRUE, ... ) 
           {
               ## Call 'plotHist()' from 'mcmcoutputfixhier'
               callNextMethod( x, dev, ... )
           }
)

setMethod( "plotDens", signature( x   = "mcmcoutputfixhierpost",
                                  dev = "ANY" ),
           function( x, dev = TRUE, ... ) 
           {
               ## Call 'plotDens()' from 'mcmcoutputfixhier'
               callNextMethod( x, dev, ... )
           }
)

setMethod("plotPointProc", signature(x      = "mcmcoutputfixhierpost",
                                     dev    = "ANY"),
          function(x, dev = TRUE, ...) 
          {
              ## Call 'plotPointProc()' from 'mcmcoutputfixhier'
              callNextMethod(x, dev, ...)
          }
)

setMethod("plotSampRep", signature(x    = "mcmcoutputfixhierpost",
                                   dev  = "ANY"),
          function(x, dev = TRUE, ...)
          {
              ## Call 'plotSampRep()' from 'mcmcoutputfixhier'
              callNextMethod(x, dev, ...)
          }
)

setMethod("plotPostDens", signature(x   = "mcmcoutputfixhierpost",
                                    dev = "ANY"),
          function(x, dev = TRUE, ...)
          {
              ## Call 'plotPostDens' from 'mcmcoutputfixhier'
              callNextMethod(x, dev, ...)
          }
)

setMethod( "subseq", signature( object    = "mcmcoutputfixhierpost",
                                index     = "array" ), 
           function( object, index ) 
           {
               ## TODO: Check arguments via .validObject ##
               dist <- object@model@dist
               ## Call 'subseq()' from 'mcmcoutputfixhier' 
               callNextMethod( object, index )
               ## post ##
               if ( dist == "poisson" ) {
                   .subseq.Poisson.Post(object, index)
               } else if ( dist == "binomial" ) {
                   .subseq.Binomial.Mcmcoutputfixpost( object, index )
               } else if ( dist %in% c( "normal", "student" ) ) {
                   .subseq.Norstud.Mcmcoutputfixpost( object, index )
               } else if ( dist %in% c( "normult", "studmult" ) ) {
                   .subseq.Normultstud.Mcmcoutputfixpost( object, index ) 
               }
           }
)

setMethod("swapElements", signature( object  = "mcmcoutputfixhierpost", 
                                     index   = "array" ),
          function(object, index ) 
          {
              if ( object@model@K == 1 ) {
                  return( object )
              } else {
                  ## Check arguments, TODO: .validObject ##
                  dist <- object@model@dist
                  ## Call 'swapElements()' from 'mcmcoutputfixhier' 
                  object <- callNextMethod( object, index )    
                  if ( dist == "poisson" ) {
                      .swapElements.Poisson.Post( object, index )
                  } else if ( dist == "binomial" ) {
                      .swapElements.Binomial.Mcmcoutputfixpost( object, index )
                  } else if ( dist %in% c( "normal", "student" ) ) {
                      .swapElements.Norstud.Mcmcoutputfixpost( object, index )
                  } else if ( dist %in% c( "normult", "studmult" ) ) {
                      .swapElements.Normultstud.Mcmcoutputfixpost( object, index ) 
                  }
              }
          }          
)
