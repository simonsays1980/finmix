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

.mcmcoutputpermfixhierpost <- setClass( "mcmcoutputpermfixhierpost",
                                        contains = c( "mcmcpermfixpost", 
                                                      "mcmcoutputfixhierpost" ),
                                        validity = function( object ) 
                                        {
                                            ## else: OK
                                            TRUE
                                        }
)

setMethod("initialize", "mcmcoutputpermfixhierpost",
          function(.Object, mcmcoutput, Mperm = integer(), 
                   parperm = list(), logperm = list(), 
                   postperm = list()) 
          {
              .Object@M         <- mcmcoutput@M
              .Object@burnin    <- mcmcoutput@burnin
              .Object@ranperm   <- mcmcoutput@ranperm
              .Object@par       <- mcmcoutput@par
              .Object@log       <- mcmcoutput@log
              .Object@hyper     <- mcmcoutput@hyper
              .Object@post      <- mcmcoutput@post
              .Object@model     <- mcmcoutput@model
              .Object@prior     <- mcmcoutput@prior
              .Object@Mperm     <- Mperm
              .Object@parperm   <- parperm
              .Object@logperm   <- logperm
              .Object@postperm  <- postperm
              .Object
          }
)

setMethod("show", "mcmcoutputpermfixhierpost",
          function(object) 
          {
              cat("Object 'mcmcoutputperm'\n")
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
              cat("     Mperm       :", object@Mperm, "\n")
              cat("     parperm     : List of", 
                  length(object@parperm), "\n")
              cat("     logperm     : List of",
                  length(object@logperm), "\n")
              cat("     postperm    : List of",
                  length(object@postperm), "\n")
              cat("     model       : Object of class", 
                  class(object@model), "\n")
              cat("     prior       : Object of class",
                  class(object@prior), "\n")
          }
)

setMethod( "plotTraces", signature( x     = "mcmcoutputpermfixhierpost", 
                                    dev   = "ANY",
                                    lik   = "ANY",
                                    col   = "ANY" ), 
          function( x, dev = TRUE, lik = 1, col = FALSE, ... ) 
          {
              dist <- x@model@dist
              if ( lik %in% c( 0, 1 ) ) {
                  if ( dist == "poisson" ) {
                      .permtraces.Poisson.Hier( x, dev )
                  } else if ( dist == "binomial" ) {
                      .permtraces.Binomial( x, dev )
                  } else if ( dist == "exponential" ) {
                      .permtraces.Exponential( x, dev )
                  } else if ( dist == "normal" ) {
                      .permtraces.Normal.Hier( x, dev )
                  } else if ( dist == "student" ) {
                      .permtraces.Student.Hier( x, dev )
                  } else if ( dist == "normult" ) {
                      .permtraces.Normult.Hier( x, dev, col )
                  } else if ( dist == "studmult" ) {
                      .permtraces.Studmult.Hier( x, dev, col )
                  }
              }    
              if ( lik %in% c( 1, 2 ) ) {
                  ## log ##
                  .permtraces.Log( x, dev )
              }
          }
)

setMethod( "plotHist", signature( x   = "mcmcoutputpermfixhierpost", 
                                  dev = "ANY" ), 
           function( x, dev = TRUE, ... ) 
           {
               dist <- x@model@dist
               if ( dist == "poisson" ) {
                   .permhist.Poisson.Hier( x, dev )
               } else if ( dist == "binomial" ) {
                   .permhist.Binomial( x, dev )
               } else if ( dist == "exponential" ) {
                   .permhist.Exponential( x, dev )
               } else if ( dist == "normal" ) {
                   .permhist.Normal.Hier( x, dev )
               } else if ( dist == "student" ) {
                   .permhist.Student.Hier( x, dev )
               } else if ( dist == "normult" ) {
                   .permhist.Normult.Hier( x, dev )
               } else if ( dist == "studmult" ) {
                   .permhist.Studmult.Hier( x, dev )
               }

           }
)

setMethod( "plotDens", signature( x   = "mcmcoutputpermfixhierpost", 
                                  dev = "ANY" ), 
           function( x, dev = TRUE, ... ) 
           {
               dist <- x@model@dist
               if ( dist == "poisson" ) {
                   .permdens.Poisson.Hier( x, dev )
               } else if ( dist == "binomial" ) {
                   .permdens.Binomial( x, dev )
               } else if ( dist == "exponential" ) {
                   .permdens.Exponential( x, dev )
               } else if ( dist == "normal" ) {
                   .permdens.Normal.Hier( x, dev )
               } else if ( dist == "student" ) {
                   .permdens.Student.Hier( x, dev )
               } else if ( dist == "normult" ) {
                   .permdens.Normult.Hier( x, dev )
               } else if ( dist == "studmult" ) {
                   .permdens.Studmult.Hier( x, dev )
               }

           }
)

setMethod( "plotPointProc", signature( x      = "mcmcoutputpermfixhierpost",
                                       dev    = "ANY" ),
            function( x, dev = TRUE, ... )
            {
                dist <- x@model@dist
                if ( dist %in% c( "poisson", "exponential" ) ) {
                    .permpointproc.Poisson( x, dev )
                } else if ( dist == "binomial" ) {
                    .permpointproc.Binomial( x, dev )
                } else if ( dist == "exponential" ) {
                    .permpointproc.Exponential( x, dev )
                } else if ( dist %in% c( "normal", "student" ) ) {
                    .permpointproc.Normal( x, dev )
                } else if ( dist %in% c( "normult", "studmult" ) ) {
                    .permpointproc.Normult( x, dev )
                }
            }
)

setMethod( "plotSampRep", signature( x    = "mcmcoutputpermfixhierpost",
                                     dev  = "ANY" ),
           function( x, dev, ... ) 
           {
               dist <- x@model@dist
               if ( dist == "poisson" ) {
                   .permsamprep.Poisson( x, dev )
               } else if ( dist == "binomial" ) {
                   .permsamprep.Binomial( x, dev )                                      
               } else if ( dist == "exponential" ) {
                   .permsamprep.Exponential( x, dev )
               } else if ( dist == "normal" ) {
                   .permsamprep.Normal( x, dev )
               } else if ( dist == "student" ) {
                   .permsamprep( x, dev )
               } else if ( dist == "normult" ) {
                   .permsamprep.Normal( x, dev)
               } else if ( dist == "studmult" ) {
                   .permsamprep.Studmult( x, dev )
               }
           }
)

setMethod( "plotPostDens", signature( x   = "mcmcoutputpermfixhierpost",
                                      dev = "ANY" ),
           function( x, dev = TRUE, ... ) 
           {
               dist <- x@model@dist
               if ( dist %in% c( "poisson", "exponential" ) ) {
                   .permpostdens.Poisson( x, dev )
               } else if ( dist == "binomial" ) {
                   .permpostdens.Binomial( x, dev )
               } else if ( dist == "normal" ) {
                   .permpostdens.Normal( x, dev )
               } else if ( dist == "student" ) {
                   .permpostdens.Student( x, dev )
               } else if ( dist == "normult" ){
                   .permpostdens.Normult( x, dev )
               } else if ( dist == "studmult" ) {
                   .permpostdens.Studmult( x, dev )
               }
           }
)

