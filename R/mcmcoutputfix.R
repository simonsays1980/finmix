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

.mcmcoutputfix <- setClass("mcmcoutputfix",
                           representation(M 	    = "integer",
                                          burnin    = "integer",
                                          ranperm   = "logical",
                                          par 	    = "list",
                                          log	    = "list",
                                          model 	= "model",
                                          prior	    = "prior"),
                           validity = function(object) 
                           {
                               ##else: OK
                               TRUE
                           },
                           prototype(M          = integer(),
                                     burnin     = integer(),
                                     ranperm    = logical(),
                                     par        = list(),
                                     log        = list(),
                                     model      = model(),
                                     prior      = prior()
                                     )
)

setMethod("show", "mcmcoutputfix", 
          function(object) {
              cat("Object 'mcmcoutputfix'\n")
              cat("     class       :", class(object), "\n")
              cat("     M           :", object@M, "\n")
              cat("     burnin      :", object@burnin, "\n")
              cat("     ranperm     :", object@ranperm, "\n")
              cat("     par         : List of", 
                  length(object@par), "\n")
              cat("     log         : List of", 
                  length(object@log), "\n")
              cat("     model       : Object of class", 
                  class(object@model), "\n")
              cat("     prior       : Object of class",
                  class(object@prior), "\n")
          }
)

setMethod( "plotTraces", signature(x     = "mcmcoutputfix", 
                                   dev   = "ANY",
                                   lik   = "ANY",
                                   col   = "ANY" ), 
           function( x, dev = TRUE, lik = 1, col = FALSE, ... ) 
           {
               dist <- x@model@dist
               if ( lik %in% c( 0, 1 ) ) {
                   if( dist == "poisson" ) {
                       .traces.Poisson( x, dev )
                   } else if ( dist == "binomial" ) {
                       .traces.Binomial( x, dev )    
                   } else if ( dist == "exponential" ) {
                      .traces.Exponential( x, dev )                      
                   } else if ( dist == "normal" ) {
                      .traces.Normal( x , dev )
                   } else if ( dist == "student" ) {
                       .traces.Student( x, dev )
                   } else if ( dist == "normult" ) {
                       .traces.Normult( x, dev, col )
                   } else if ( dist == "studmult" ) {
                       .traces.Studmult( x, dev, col )
                   }
               }
               if ( lik %in% c( 1, 2 ) ) {
                  ## log ##
                  .traces.Log( x, dev, col )
               }
           }
)

setMethod( "plotHist", signature( x   = "mcmcoutputfix", 
                                  dev = "ANY" ), 
          function( x, dev = TRUE, ... ) 
          {
              dist <- x@model@dist
              if( dist == "poisson" ) {
                  .hist.Poisson( x, dev )
              } else if ( dist == "binomial" ) {
                  .hist.Binomial( x, dev )
              } else if ( dist == "exponential" ) {
                  .hist.Exponential( x, dev )
              } else if ( dist == "normal" ) {
                  .hist.Normal( x, dev )
              } else if ( dist == "student" ) {
                  .hist.Student( x, dev )
              } else if ( dist == "normult" ) {
                  .hist.Normult( x, dev )
              } else if ( dist == "studmult" ) {
                  .hist.Studmult( x, dev )
              }
          }
)

setMethod( "plotDens", signature( x   = "mcmcoutputfix",
                                  dev = "ANY" ),
          function( x, dev = TRUE, ... )
          {
              dist <- x@model@dist
              if ( dist == "poisson" ) {
                  .dens.Poisson( x, dev )
              } else if ( dist == "binomial" ) {
                  .dens.Binomial( x, dev )
              } else if ( dist == "exponential" ) {
                  .dens.Exponential( x, dev )
              } else if ( dist == "normal" ) {
                  .dens.Normal( x, dev )
              } else if ( dist == "student" ) {
                  .dens.Student( x, dev )
              } else if ( dist == "normult" ) {
                  .dens.Normult( x, dev )
              } else if ( dist == "studmult" ) {
                  .dens.Studmult( x, dev )
              }
          }
)

setMethod("plotPointProc", signature( x      = "mcmcoutputfix",
                                      dev    = "ANY"),
          function(x, dev = TRUE, ...)
          {
              dist <- x@model@dist
              if (dist == "poisson") {
                  .pointproc.Poisson(x, dev)
              } else if (dist == "binomial") {
                  .pointproc.Binomial(x, dev)
              } 
          }
)

setMethod("plotSampRep", signature(x    = "mcmcoutputfix",
                                   dev  = "ANY"),
          function(x, dev, ...) 
          {
              dist <- x@model@dist
              if (dist == "poisson") {
                  .samprep.Poisson(x, dev)
              } else if (dist == "binomial") {
                  .samprep.Binomial(x, dev)
              }
          }
)

setMethod("plotPostDens", signature(x   = "mcmcoutputfix",
                                    dev = "ANY"),
          function(x, dev = TRUE, ...) 
          {
              dist <- x@model@dist
              if (dist == "poisson") {
                  .postdens.Poisson(x, dev)
              } else if (dist == "binomial") {
                  .postdens.Binomial(x, dev)
              }
          }
)

setMethod( "subseq", signature( object    = "mcmcoutputfix", 
                                index     = "array" ), 
          function( object, index ) 
          {
              .subseq.valid.Arg( object, index )
              dist      <- object@model@dist
              object@M  <- sum( index )
              ## log ##
              object    <- .subseq.Log.Fix( object, index )
              ## par ##
              if ( dist == "poisson" ) {
                  .subseq.Poisson( object, index )
              } else if ( dist == "binomial" ) {
                  .subseq.Binomial( object, index )
              } else if ( dist == "exponential" ) {
                  .subseq.Exponential ( object, index )              
              } else if ( dist == "normal" ) {
                  .subseq.Normal( object, index )
              } else if ( dist == "student" ) {
                  .subseq.Student( object, index )
              } else if ( dist == "normult" ) {
                  .subseq.Normult( object, index )
              } else if ( dist == "studmult" ) {
                  .subseq.Studmult( object, index )
              }
          }
)

setMethod( "swapElements", signature( object  = "mcmcoutputfix", 
                                      index   = "array" ),
          function( object, index ) 
          { ## Check arguments, TODO: .validObject ##
              .swapElements.valid.Arg( object, index )
              if ( object@model@K == 1 ) {
                  return( object )
              } else {
                  dist <- object@model@dist
                  if ( dist == "poisson" ) {
                      .swapElements.Poisson( object, index )
                  } else if ( dist == "binomial" ) {
                      .swapElements.Binomial( object, index )
                  } else if ( dist == "exponential" ) {
                      .swapElements.Exponential( object, index )
                  } else if ( dist == "normal" ) {
                      .swapElements.Normal( object, index )
                  } else if ( dist == "student" ) {
                      .swapElements.Student( object, index )
                  } else if ( dist == "normult" ) {
                      .swapElements.Normult( object, index )
                  } else if ( dist == "studmult" ) {
                      .swapElements.Studmult( object, index )
                  }
              }
          }
) 
          
setMethod( "extract", signature( object = "mcmcoutputfix",
                                 index  = "numeric" ),
           function( object, index ) 
           {
               dist <- object@model@dist
               if ( dist == "normult" ) {
                   .extract.Normult( object, as.integer( index ) ) 
               }
           }
)

setMethod( "moments", signature( object = "mcmcoutputfix" ),
           function( object ) 
           {
               dist <- objject@model@dist
               if ( dist == "normult" ) {
                   .moments.Normult.Mcmcoutput( object )
               }
           }
)
## Getters ##
setMethod( "getM", "mcmcoutputfix", 
          function( object ) 
          {
              return( object@M )
          }
)

setMethod("getBurnin", "mcmcoutputfix",
          function(object) 
          {
              return(object@burnin)
          }
)

setMethod("getRanperm", "mcmcoutputfix", 
          function(object) 
          {
              return(object@ranperm)
          }
)

setMethod("getPar", "mcmcoutputfix", 
          function(object) 
          {
              return(object@par)
          }
)

setMethod("getLog", "mcmcoutputfix", 
          function(object) 
          {
              return(object@log)
          }
)

setMethod("getModel", "mcmcoutputfix", 
          function(object) 
          {
              return(object@model)
          }
)

setMethod("getPrior", "mcmcoutputfix", 
          function(object) 
          {
              return(object@prior)
          }
)

## No setters as users are not intended to manipulate ##
## this object ##i

### Private functions
### These functions are not exported

### Plot
### Traces Poisson: Plots the traces of MCMC samples
### for Poisson mixture. If dev = FALSE, no graphical
### device is started, instead it is assumed that the
### user wants to save the graphic to a file.
".traces.Poisson" <- function(x, dev)
{
    K <- x@model@K
    trace.n <- K
    if (.check.grDevice() && dev) {	
        dev.new(title = "Traceplots")
    }
    par(mfrow = c(trace.n, 1), mar = c(1, 0, 0, 0), 
        oma = c(4, 5, 4, 4))
    lambda <- x@par$lambda
    for (k in 1:K) {
        plot(lambda[, k], type = "l", axes = F, 
             col = "gray20", xlab = "", ylab = "")
        axis(2, las = 2, cex.axis = .7)
        mtext(side = 2, las = 2, bquote(lambda[k = .(k)]), 
              cex = .6, line = 3)
    }
    axis(1)
    mtext(side = 1, "Iterations", cex = .7, line = 3)		
}

".traces.Binomial" <- function(x, dev)
{
    K <- x@model@K
    trace.n <- K
    if (.check.grDevice() && dev) {
        dev.new(title = "Traceplots")        
    } 
    par(mfrow = c(trace.n, 1), mar = c(1, 0, 0, 0),
        oma = c(4, 5, 4, 4))
    p <- x@par$p
    for (k in 1:K) {
        plot(p[, k], type = "l", axes = F, col = "gray20",
             xlab = "", ylab = "")
        axis(2, las = 2, cex.axis = .7)
        mtext(side = 2, las = 2, bquote(p[k = .(k)]),
              cex = .6, line = 3)
    }
    axis(1)
    mtext(side = 1, "Iterations", cex = .7, line = 3)
}

### --------------------------------------------------------------------
### .traces.Exponential
### @description    Plots traces for parameters of Exponential mixture.
### @par    x       an object of class mcmcoutputfix
###         dev     an object of class 'logical' 
### @detail         Plots the traces for each component parameter of an
###                 Exponential mixture. If 'dev' is set to FALSE 
###                 (TRUE is default) no device is created, instead 
###                 the graphic can be stored to a file.
### @see            ?mcmcoutput, ?plotTraces
### @author         Lars Simon Zehnder
### --------------------------------------------------------------------
".traces.Exponential" <- function( x, dev )
{
    K       <- x@model@K
    trace.n <- K
    if ( .check.grDevice() && dev ) {	
        dev.new( title = "Traceplots" )
    }
    par( mfrow = c( trace.n, 1 ), mar = c( 1, 0, 0, 0 ), 
         oma = c( 4, 5, 4, 4 ) )
    lambda <- x@par$lambda
    for ( k in 1:K ) {
        plot( lambda[, k], type = "l", axes = F, 
              col = "gray20", xlab = "", ylab = "" )
        axis( 2, las = 2, cex.axis = .7 )
        mtext( side = 2, las = 2, bquote( lambda[k = .( k )] ), 
               cex = .6, line = 3 )
    }
    axis( 1 )
    mtext( side = 1, "Iterations", cex = .7, line = 3 )		
}

### --------------------------------------------------------------------
### .traces.Normal
### @description    Plots traces for parameters of a univariate Normal 
###                 mixture.
### @par    x       an object of class mcmcoutputfix
###         dev     an object of class 'logical' 
### @detail         Plots the traces for each component parameter of an
###                 Normal mixture. If 'dev' is set to FALSE 
###                 (TRUE is default) no device is created, instead 
###                 the graphic can be stored to a file.
### @see            ?mcmcoutput, ?plotTraces
### @author         Lars Simon Zehnder
### --------------------------------------------------------------------
".traces.Normal" <- function( x, dev ) 
{
    K       <- x@model@K
    trace.n <- 2 * K
    if ( .check.grDevice() && dev ) {
        dev.new( title = "Traceplots" )        
    }
    par( mfrow = c( trace.n, 1 ), mar = c( 1, 0, 0, 0 ),
         oma = c( 4, 5, 4, 4 ) )
    mu      <- x@par$mu
    sigma   <- x@par$sigma
    for ( k in 1:K ) {
        plot( mu[, k], type = "l", axes = F,
             col = "gray20", xlab = "", ylab = "" )
        axis( 2, las = 2, cex.axis = .7 )
        mtext( side = 2, las = 2, bquote( mu[k = .( k )] ),
               cex = .6, line = 3 )       
    }
    for ( k in 1:K ) {
        plot( sigma[, k], type = "l", axes = F,
              col = "gray30", xlab = "", ylab = "" )
        axis( 2, las = 2, cex.axis = .7 )              
        mtext( side = 2, las = 2, bquote( sigma[k = .( k )]),
               cex = .6, line = 3 )        
    }
    axis( 1 ) 
    mtext( side = 1, "Iterations", cex = .7, line = 3 )
}

### --------------------------------------------------------------------
### .traces.Student
### @description    Plots traces for parameters of a univariate Student 
###                 mixture.
### @par    x       an object of class mcmcoutputfix
###         dev     an object of class 'logical' 
### @detail         Plots the traces for each component parameter of an
###                 Student mixture. If 'dev' is set to FALSE 
###                 (TRUE is default) no device is created, instead 
###                 the graphic can be stored to a file.
### @see            ?mcmcoutput, ?plotTraces
### @author         Lars Simon Zehnder
### --------------------------------------------------------------------
".traces.Student" <- function( x, dev ) 
{
    K       <- x@model@K
    trace.n <- 3 * K
    if ( .check.grDevice() && dev ) {
        dev.new( title = "Traceplots" )        
    }
    par( mfrow = c( trace.n, 1 ), mar = c( 1, 0, 0, 0 ),
         oma = c( 4, 5, 4, 4 ) )
    mu      <- x@par$mu
    sigma   <- x@par$sigma
    df      <- x@par$df
    for ( k in 1:K ) {
        plot( mu[, k], type = "l", axes = F,
             col = "gray20", xlab = "", ylab = "" )
        axis( 2, las = 2, cex.axis = .7 )
        mtext( side = 2, las = 2, bquote( mu[k = .( k )] ),
               cex = .6, line = 3 )       
    }
    for ( k in 1:K ) {
        plot( sigma[, k], type = "l", axes = F,
              col = "gray30", xlab = "", ylab = "" )
        axis( 2, las = 2, cex.axis = .7 )              
        mtext( side = 2, las = 2, bquote( sigma[k = .( k )]),
               cex = .6, line = 3 )        
    }
    for ( k in 1:K ) {
        plot( df[, k], type = "l", axes = F,
              col = "gray40", xlab = "", ylab = "" )
        axis( 2, las = 2, cex.axis = .7 )
        mtext( side = 2, las = 2, bquote( nu[k = .( k )]),
               cex = .6, line = 3 )
    }
    axis( 1 ) 
    mtext( side = 1, "Iterations", cex = .7, line = 3 )
}

### -------------------------------------------------------------------------
### .traces.Normult
### @description    Plots the traces of parameters and moments of a multi-
###                 variate Normal distribution.
### @par    x       an mcmcoutputfix object
###         dev     a logical
###         col     a logical
### @return         a graphical device
### @detail         If dev = FALSE, the plot can be sent to a file. In case
###                 col = TRUE, rainbow colors are used.
### @see            ?plotTraces
### @author Lars Simon Zehnder
### -------------------------------------------------------------------------
".traces.Normult" <- function( x, dev, col ) 
{
    K       <- x@model@K
    r       <- x@model@r
    if ( .check.grDevice() && dev ) {
        dev.new( title = "Traceplots" ) 
    }
    trace.n <- r + 2
    par( mfrow = c( trace.n, 1 ), mar = c( 1, 2, 0, 0 ),
         oma = c( 4, 5, 2, 4 ) )
    mu      <- x@par$mu
    sigma   <- x@par$sigma
    if ( col ) {
        cscale  <- rainbow( K, start = 0.5, end = 0 )
    } else {
        cscale  <- gray.colors( K, start = 0.5, end = 0.15 )
    }
    for ( rr in 1:r ) {
        mmax    <- max( mu[,rr,] )
        mmin    <- min( mu[,rr,] )
        plot( mu[, rr, 1], type = "l", axes = F, 
             col = cscale[1], xlab = "", ylab = "",
             ylim = c( mmin, mmax + 0.3 * (mmax - mmin) ) )
        for ( k in 2:K ) {
            lines( mu[, rr, k], col = cscale[ k ] )
        }
        axis( 2, las = 2, cex.axis = .7 )
        mtext( side = 2, las = 2, bquote( mu[rr = .( rr )] ),
               cex = .6, line = 3 )       
        if ( rr == 1 ) {
            name    <- vector( "character", K )
            for (k in 1:K ) {
                name[k] <- paste( "k = ", k, sep = "")
            }
            legend( "top", legend = name, col = cscale, horiz = TRUE, 
                    lty = 1 )
        }
    }
    sigma.tr    <- array( numeric(), dim = c( x@M, K ) )
    sigma.det   <- array( numeric(), dim = c( x@M, K ) )                         
    for ( k in 1:K ) {
        sigma.tr[, k]    <- sapply( seq( 1, x@M ), 
                                    function( i ) sum( diag( qinmatr( sigma[i,, k] ) ) ) )
        sigma.det[, k]   <- sapply( seq( 1, x@M ),
                                    function( i ) log( det( qinmatr( sigma[i,, k] ) ) ) )        
    }
    # Sigma traces 
    mmax    <- max( sigma.tr ) 
    mmin    <- min( sigma.tr )
    plot( sigma.tr[, 1], type = "l", axes = F, 
          col = cscale[1], xlab = "", ylab = "",
          ylim = c( mmin, mmax ) )
    for ( k in 2:K ) {
        lines( sigma.tr[, k], col = cscale[k] )
    }    
    axis( 2, las = 2, cex.axis = .7 )
    mtext( side = 2, las = 2, bquote( tr( Sigma ) ),
           cex = .6, line = 3 )        
    
    # Sigma logdets
    mmax    <- max( sigma.det )
    mmin    <- min( sigma.det )
    plot( sigma.det[, 1], type = "l", axes = F,
          col = cscale[1], xlab = "", ylab = "",
          ylim = c( mmin, mmax ) )
    for( k in 2:K ) {
        lines( sigma.det[, k], col = cscale[k] )
    }
    axis( 2, las = 2, cex.axis = .7 )
    mtext( side = 2, las = 2, bquote( log( det( Sigma ) ) ),
           cex = .6, line = 3 )

    axis( 1 ) 
    mtext( side = 1, "Iterations", cex = .7, line = 3 )

    # Get moments 
    moms    <- moments_cc( x )
    for ( rr in 1:r ) {
        if ( .check.grDevice() && dev ) {
            dev.new( title = paste( "Traceplots Feature ", rr, sep = "" ) )
        }
        par( mfrow = c( 2, 2 ), mar = c( 4, 4, 0.5, 0.5 ), 
             oma = c( 1.5, 2, 1, 1 ) )
        # Mu
        plot( moms$mean[, rr], type = "l", axes = F,
             xlab = "", ylab = "", col = cscale[K] )
        axis( 2, las = 2, cex.axis = .7 )
        mtext( side = 2, las = 2, bquote( mu ), cex = .6,
               line = 3 )
        axis( 1 )
        mtext( side = 1, "Iterations", cex = .7, line = 3 )       
        # Variance
        plot( moms$var[, rr], type = "l", axes = F,
             xlab = "", ylab = "", col = cscale[K] )
        axis( 2, las = 2, cex.axis = .7 )
        mtext( side = 2, las = 2, bquote( sigma ), cex = .6,
               line = 3 )
        axis( 1 )
        mtext( side = 1, "Iterations", cex = .7, line = 3 )       
        # Skewness
        plot( moms$skewness[, rr], type = "l", axes = F,
             xlab = "", ylab = "", col = cscale[K] )
        axis( 2, las = 2, cex.axis = .7 )
        mtext( side = 2, "Skewness", cex = .6,
               line = 3 )
        axis( 1 )
        mtext( side = 1, "Iterations", cex = .7, line = 3 )       
        # Kurtosis
        plot( moms$kurtosis[, rr], type = "l", axes = F,
             xlab = "", ylab = "", col = cscale[K] )
        axis( 2, las = 2, cex.axis = .7 )
        mtext( side = 2, "Kurtosis", cex = .6,
               line = 3 )
        axis( 1 )
        mtext( side = 1, "Iterations", cex = .7, line = 3 )        
    }

}

### -------------------------------------------------------------------------
### .traces.Studmult
### @description    Plots the traces of parameters and moments of a multi-
###                 variate Student-t distribution.
### @par    x       an mcmcoutputfix object
###         dev     a logical
###         col     a logical
### @return         a graphical device
### @detail         If dev = FALSE, the plot can be sent to a file. In case
###                 col = TRUE, rainbow colors are used.
### @see            ?plotTraces
### @author Lars Simon Zehnder
### -------------------------------------------------------------------------
".traces.Studmult" <- function( x, dev, col ) 
{
    K       <- x@model@K
    r       <- x@model@r
    if ( .check.grDevice() && dev ) {
        dev.new( title = "Traceplots" ) 
    }
    trace.n <- r + 2
    par( mfrow = c( trace.n, 1 ), mar = c( 1, 2, 0, 0 ),
         oma = c( 4, 5, 4, 4 ) )
    mu      <- x@par$mu
    sigma   <- x@par$sigma
    if ( col ) {
        cscale  <- rainbow( K, start = 0.5, end = 0 )
    } else {
        cscale  <- gray.colors( K, start = 0.5, end = 0.15 )
    }
    for ( rr in 1:r ) {
        mmax    <- max( mu[,rr,] )
        mmin    <- min( mu[,rr,] )
        plot( mu[, rr, 1], type = "l", axes = F, 
             col = cscale[1], xlab = "", ylab = "",
             ylim = c( mmin, mmax + 0.3 * ( mmax  - mmin ) ) )
        for ( k in 2:K ) {
            lines( mu[, rr, k], col = cscale[ k ] )
        }
        axis( 2, las = 2, cex.axis = .7 )
        mtext( side = 2, las = 2, bquote( mu[rr = .( rr )] ),
               cex = .6, line = 3 ) 
        if ( rr == 1 ) {
            name    <- vector( "character", K )
            for (k in 1:K ) {
                name[k] <- paste( "k = ", k, sep = "")
            }
            legend( "top", legend = name, col = cscale, horiz = TRUE, 
                    lty = 1 )
        }
    }
    sigma.tr    <- array( numeric(), dim = c( x@M, K ) )
    sigma.det   <- array( numeric(), dim = c( x@M, K ) )                         
    for ( k in 1:K ) {
        sigma.tr[, k]    <- sapply( seq( 1, x@M ), 
                                    function( i ) sum( diag( qinmatr( sigma[i,, k] ) ) ) )
        sigma.det[, k]   <- sapply( seq( 1, x@M ),
                                    function( i ) log( det( qinmatr( sigma[i,, k] ) ) ) )        
    }
    # Sigma traces 
    mmax    <- max( sigma.tr ) 
    mmin    <- min( sigma.tr )
    plot( sigma.tr[, 1], type = "l", axes = F, 
          col = cscale[1], xlab = "", ylab = "",
          ylim = c( mmin, mmax ) )
    for ( k in 2:K ) {
        lines( sigma.tr[, k], col = cscale[k] )
    }
    axis( 2, las = 2, cex.axis = .7 )
    mtext( side = 2, las = 2, bquote( tr( Sigma ) ),
           cex = .6, line = 3 )        
    
    # Sigma logdets
    mmax    <- max( sigma.det )
    mmin    <- min( sigma.det )
    plot( sigma.det[, 1], type = "l", axes = F,
          col = cscale[1], xlab = "", ylab = "",
          ylim = c( mmin, mmax ) )
    for( k in 2:K ) {
        lines( sigma.det[, k], col = cscale[k] )
    }
    axis( 2, las = 2, cex.axis = .7 )
    mtext( side = 2, las = 2, bquote( log( det( Sigma ) ) ),
           cex = .6, line = 3 )

    axis( 1 ) 
    mtext( side = 1, "Iterations", cex = .7, line = 3 )
    
    # Degrees of freedom
    if ( .check.grDevice() && dev ) {
        dev.new( title = "Traceplots Degrees of Freedom" )
    }
    degf    <- x@par$df  
    par( mfrow = c( K, 1 ), mar = c( 1, 2, 0, 0),
         oma = c( 4, 5, 4, 4 ) )
    for ( k in 1:K ) {        
        plot( degf[,k], type = "l", axes = F,
              col = cscale[K], xlab = "", ylab = "" )
        axis( 2, las = 2, cex.axis = .7 )
        mtext( side = 2, las = 2, bquote( nu[ k = .(k) ] ),
               cex = .6, line = 3 )
    } 
    axis( 1 ) 
    mtext( side = 1, "Iterations", cex = .7, line = 3 )
    # Get moments 
    moms    <- moments_cc( x )
    for ( rr in 1:r ) {
        if ( .check.grDevice() && dev ) {
            dev.new( title = paste( "Traceplots Feature ", rr, sep = "" ) )
        }
        par( mfrow = c( 2, 2 ), mar = c( 4, 4, 0.5, 0.5 ), 
             oma = c( 1.5, 2, 1, 1 ) )
        # Mu
        plot( moms$mean[, rr], type = "l", axes = F,
             xlab = "", ylab = "", col = cscale[K] )
        axis( 2, las = 2, cex.axis = .7 )
        mtext( side = 2, las = 2, bquote( mu ), cex = .6,
               line = 3 )
        axis( 1 )
        mtext( side = 1, "Iterations", cex = .7, line = 3 )       
        # Variance
        plot( moms$var[, rr], type = "l", axes = F,
             xlab = "", ylab = "", col = cscale[K] )
        axis( 2, las = 2, cex.axis = .7 )
        mtext( side = 2, las = 2, bquote( sigma ), cex = .6,
               line = 3 )
        axis( 1 )
        mtext( side = 1, "Iterations", cex = .7, line = 3 )       
        # Skewness
        plot( moms$skewness[, rr], type = "l", axes = F,
             xlab = "", ylab = "", col = cscale[K] )
        axis( 2, las = 2, cex.axis = .7 )
        mtext( side = 2, "Skewness", cex = .6,
               line = 3 )
        axis( 1 )
        mtext( side = 1, "Iterations", cex = .7, line = 3 )       
        # Kurtosis
        plot( moms$kurtosis[, rr], type = "l", axes = F,
             xlab = "", ylab = "", col = cscale[K] )
        axis( 2, las = 2, cex.axis = .7 )
        mtext( side = 2, "Kurtosis", cex = .6,
               line = 3 )
        axis( 1 )
        mtext( side = 1, "Iterations", cex = .7, line = 3 )        
    }
}

### Traces Poisson: Plots the traces of MCMC samples
### for the log-likelihoods. If dev = FALSE, no graphical
### device is started, instead it is assumed that the
### user wants to save the graphic to a file.
".traces.Log" <- function( x, dev, col )
{
    if( .check.grDevice() && dev ) {
        dev.new( title = "Log Likelihood Traceplots" )
    }
    if ( col ) {
        cscale  <- rainbow( 3, start = 0.5, end = 0 )        
    } else {
        cscale  <- gray.colors( 3, start = 0.5, end = 0.15 )
    }
    par( mfrow = c( 2, 1 ), mar = c( 1, 0, 0, 0 ),
         oma = c( 4, 5, 4, 4 ) )
    mixlik <- x@log$mixlik
    plot( mixlik, type = "l", axes = F,
          col = cscale[3], xlab = "", ylab = "" )
    axis( 2, las = 2, cex.axis = 0.7 )
    mtext( side = 2, las = 3, "mixlik", cex = 0.6,
           line = 3 )
    mixprior <- x@log$mixprior
    plot( mixprior, type = "l", axes = F,
          col = cscale[2], xlab = "", ylab = "" )
    axis( 2, las = 2, cex.axis = 0.7 )
    mtext( side = 2, las = 3, "mixprior", cex = 0.6,
           line = 3 )
    axis( 1 )
    mtext( side = 1, "Iterations", cex = 0.7, line = 3 )		
}

### Plot Histogramms
### Plot Hist Poisson: Plots Histograms for each component
### parameter. 
".hist.Poisson" <- function(x, dev)
{
    K <- x@model@K 
    if (.check.grDevice() && dev) {
        dev.new(title = "Histograms")
    }
    lambda <- x@par$lambda
    if (K == 1) {
        .symmetric.Hist(lambda, list(bquote(lambda)))
    } else {
        lab.names <- vector("list", K)
        for (k in 1:K) {
            lab.names[[k]] <- bquote(lambda[.(k)])
        }
        .symmetric.Hist(lambda, lab.names)
    }
}

".hist.Binomial" <- function(x, dev)
{
    K <- x@model@K 
    if (.check.grDevice() && dev) {
        dev.new(title = "Histograms")
    }
    p <- x@par$p
    if (K == 1) {
        .symmetric.Hist(p, list(bquote(p)))
    } else {
        lab.names <- vector("list", K)
        for (k in 1:K) {
            lab.names[[k]] <- bquote(p[.(k)])
        }
        .symmetric.Hist(p, lab.names)
    }
}

".hist.Exponential" <- function( x, dev )
{
    K <- x@model@K 
    if ( .check.grDevice() && dev ) {
        dev.new( title = "Histograms" )
    }
    lambda <- x@par$lambda
    if ( K == 1 ) {
        .symmetric.Hist( lambda, list( bquote( lambda ) ) )
    } else {
        lab.names <- vector( "list", K )
        for ( k in 1:K ) {
            lab.names[[k]] <- bquote( lambda[.( k )] )
        }
        .symmetric.Hist( lambda, lab.names )
    }
}

".hist.Normal" <- function( x, dev )
{
    K <- x@model@K 
    mu      <- x@par$mu
    sigma   <- x@par$sigma
    if ( K == 1 ) {
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Histogram Mu" )
        }
        .symmetric.Hist( mu, list( bquote( mu ) ) )
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Histogram Sigma" )
        }
        .symmetric.Hist( sigma, list( bquote( sigma ) ) ) 
    } else {
        mu.lab.names    <- vector( "list", K )
        sigma.lab.names <- vector( "list", K )
        for ( k in 1:K ) {
            mu.lab.names[[k]]       <- bquote( mu[.( k )] )
            sigma.lab.names[[k]]    <- bquote( sigma[.( k )] )
        }
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Histograms Mu" )
        }
        .symmetric.Hist( mu, mu.lab.names )
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Histograms Sigma" )
        }
        .symmetric.Hist( sigma, sigma.lab.names )
    }
}

".hist.Student" <- function( x, dev )
{
    K <- x@model@K 
    mu      <- x@par$mu
    sigma   <- x@par$sigma
    degf    <- x@par$df
    if ( K == 1 ) {
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Histogram Mu" )
        }
        .symmetric.Hist( mu, list( bquote( mu ) ) )
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Histogram Sigma" )
        }
        .symmetric.Hist( sigma, list( bquote( sigma ) ) )
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Histogram Degrees of Freedom" )
        }
        .symmetric.Hist( degf, list( bquote( nu ) ) )
    } else {
        mu.lab.names    <- vector( "list", K )
        sigma.lab.names <- vector( "list", K )
        degf.lab.names  <- vector( "list", K )
        for ( k in 1:K ) {
            mu.lab.names[[k]]       <- bquote( mu[.( k )] )
            sigma.lab.names[[k]]    <- bquote( sigma[.( k )] )
            degf.lab.names[[k]]     <- bquote( nu[.( k )] )
        }
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Histograms Mu" )
        }
        .symmetric.Hist( mu, mu.lab.names )
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Histograms Sigma" )
        }
        .symmetric.Hist( sigma, sigma.lab.names )
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Histograms Degrees of Freedom" )
        }
        .symmetric.Hist( degf, degf.lab.names )
    }
}

".hist.Normult"  <- function( x, dev ) 
{
    K       <- x@model@K
    r       <- x@model@r
    mu      <- x@par$mu
    sigma   <- x@par$sigma
    for ( rr in 1:r ) {
        if ( K == 1 ) {
            if (.check.grDevice() & dev ) {
                dev.new( title = paste( "Histograms Feature ", rr, 
                                        " Mu", sep = "" ) )
            }
            .symmetric.Hist( mu[, rr,], list( bquote( mu ) ) )
            if (.check.grDevice() & dev ) {
                dev.new( title = paste( "Histograms Feature ", rr, 
                                        " Sigma", sep = "" ) )
            }
            .symmetric.Hist( sigma[, rr,], list( bquote( sigma ) ) )           
        } else {
            mu.lab.names    <- vector( "list", K )
            sigma.lab.names <- vector( "list", K ) 
            for ( k in 1:K ) {
                mu.lab.names[[k]]       <- bquote( mu[.( k )] )
                sigma.lab.names[[k]]    <- bquote( sigma[.( k )] )
            }
            if (.check.grDevice() & dev ) {
                dev.new( title = paste( "Histograms Feature ", rr, 
                                        " Mu", sep = "" ) )
            }
            .symmetric.Hist( mu[, rr,], mu.lab.names )
            if (.check.grDevice() & dev ) {
                dev.new( title = paste( "Histograms Feature ", rr, 
                                        " Sigma", sep = "" ) )
            }
            .symmetric.Hist( sigma[, rr,], sigma.lab.names )
        }        
    }
}

".hist.Studmult"  <- function( x, dev ) 
{
    K       <- x@model@K
    r       <- x@model@r
    mu      <- x@par$mu
    sigma   <- x@par$sigma
    degf    <- x@par$df
    for ( rr in 1:r ) {
        if ( K == 1 ) {
            if (.check.grDevice() & dev ) {
                dev.new( title = paste( "Histograms Feature ", rr, 
                                        " Mu", sep = "" ) )
            }
            .symmetric.Hist( mu[, rr,], list( bquote( mu ) ) )
            if (.check.grDevice() & dev ) {
                dev.new( title = paste( "Histograms Feature ", rr, 
                                        " Sigma", sep = "" ) )
            }
            .symmetric.Hist( sigma[, rr,], list( bquote( sigma ) ) )
         
        } else {
            mu.lab.names    <- vector( "list", K )
            sigma.lab.names <- vector( "list", K ) 
            for ( k in 1:K ) {
                mu.lab.names[[k]]       <- bquote( mu[.( k )] )
                sigma.lab.names[[k]]    <- bquote( sigma[.( k )] )
            }
            if (.check.grDevice() & dev ) {
                dev.new( title = paste( "Histograms Feature ", rr, 
                                        " Mu", sep = "" ) )
            }
            .symmetric.Hist( mu[, rr,], mu.lab.names )
            if (.check.grDevice() & dev ) {
                dev.new( title = paste( "Histograms Feature ", rr, 
                                        " Sigma", sep = "" ) )
            }
            .symmetric.Hist( sigma[, rr,], sigma.lab.names )
        }
    }
    if ( K == 1 ) { 
        if (.check.grDevice() & dev ) {
            dev.new( title = paste( "Histograms Feature ", rr, 
                                    " Mu", sep = "" ) )
        }
        .symmetric.Hist( degf[, rr,], list( bquote( nu ) ) ) 
    } else {
        degf.lab.names  <- vector( "list", K )
        for ( k in 1:K ) {
            degf.lab.names[[k]] <- bquote( nu[.( k )] )
        }
        if (.check.grDevice() & dev ) {
            dev.new( title = paste( "Histograms Feature ", rr, 
                                    " Sigma", sep = "" ) )
        }
        .symmetric.Hist( degf[, rr,], degf.lab.names ) 
    }
}

### Plot Densities
### Plot Dens Poisson: Plots Kernel densities for each
### component parameter.
".dens.Poisson" <- function(x, dev)
{
    K   <- x@model@K
    if (.check.grDevice() && dev) {
        dev.new(title = "Densities")
    }
    lambda  <- x@par$lambda
    if (K == 1) {
        .symmetric.Dens(lambda, list(bquote(lambda)))
    } else {
        lab.names   <- vector("list", K)
        for (k in seq(1, K)) {
            lab.names[[k]]  <- bquote(lambda[.(k)])
        }
        .symmetric.Dens(lambda, lab.names)
    }
}

".dens.Binomial" <- function(x, dev)
{
    K   <- x@model@K
    if (.check.grDevice() && dev) {
        dev.new(title = "Densities")
    }
    p  <- x@par$p
    if (K == 1) {
        .symmetric.Dens(p, list(bquote(p)))
    } else {
        lab.names   <- vector("list", K)
        for (k in seq(1, K)) {
            lab.names[[k]]  <- bquote(p[.(k)])
        }
        .symmetric.Dens(p, lab.names)
    }
}

".dens.Exponential" <- function( x, dev )
{
    K <- x@model@K 
    if ( .check.grDevice() && dev ) {
        dev.new( title = "Densities" )
    }
    lambda <- x@par$lambda
    if ( K == 1 ) {
        .symmetric.Dens( lambda, list( bquote( lambda ) ) )
    } else {
        lab.names <- vector( "list", K )
        for ( k in 1:K ) {
            lab.names[[k]] <- bquote( lambda[.( k )] )
        }
        .symmetric.Dens( lambda, lab.names )
    }
}

".dens.Normal" <- function( x, dev )
{
    K <- x@model@K 
    mu      <- x@par$mu
    sigma   <- x@par$sigma
    if ( K == 1 ) {
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Density Mu" )
        }
        .symmetric.Dens( mu, list( bquote( mu ) ) )
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Density Sigma" )
        }
        .symmetric.Dens( sigma, list( bquote( sigma ) ) ) 
    } else {
        mu.lab.names    <- vector( "list", K )
        sigma.lab.names <- vector( "list", K )
        for ( k in 1:K ) {
            mu.lab.names[[k]]       <- bquote( mu[.( k )] )
            sigma.lab.names[[k]]    <- bquote( sigma[.( k )] )
        }
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Densities Mu" )
        }
        .symmetric.Dens( mu, mu.lab.names )
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Densities Sigma" )
        }
        .symmetric.Dens( sigma, sigma.lab.names )
    }
}

".dens.Student.Hier" <- function( x, dev )
{
    K <- x@model@K 
    mu      <- x@par$mu
    sigma   <- x@par$sigma
    degf    <- x@par$df
    if ( K == 1 ) {
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Density Mu" )
        }
        .symmetric.Dens( mu, list( bquote( mu ) ) )
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Density Sigma" )
        }
        .symmetric.Dens( sigma, list( bquote( sigma ) ) )
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Density Degrees of Freedom" )
        }
        .symmetric.Dens( degf, list( bquote( nu ) ) )
    } else {
        mu.lab.names    <- vector( "list", K )
        sigma.lab.names <- vector( "list", K )
        degf.lab.names  <- vector( "list", K )
        for ( k in 1:K ) {
            mu.lab.names[[k]]       <- bquote( mu[.( k )] )
            sigma.lab.names[[k]]    <- bquote( sigma[.( k )] )
            degf.lab.names[[k]]     <- bquote( nu[.( k )] )
        }
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Densities Mu" )
        }
        .symmetric.Dens( mu, mu.lab.names )
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Densities Sigma" )
        }
        .symmetric.Dens( sigma, sigma.lab.names )
        if ( .check.grDevice() && dev ) {
            dev.new( title = "Densities Degrees of Freedom" )
        }
        .symmetric.Dens( degf, degf.lab.names )
    }
}

".dens.Normult"  <- function( x, dev ) 
{
    K       <- x@model@K
    r       <- x@model@r
    mu      <- x@par$mu
    sigma   <- x@par$sigma
    for ( rr in 1:r ) {
        if ( K == 1 ) {
            if (.check.grDevice() & dev ) {
                dev.new( title = paste( "Densities Feature ", rr, 
                                        " Mu", sep = "" ) )
            }
            .symmetric.Dens( mu[, rr,], list( bquote( mu ) ) )
            if (.check.grDevice() & dev ) {
                dev.new( title = paste( "Densities Feature ", rr, 
                                        " Sigma", sep = "" ) )
            }
            .symmetric.Dens( sigma[, rr,], list( bquote( sigma ) ) )           
        } else {
            mu.lab.names    <- vector( "list", K )
            sigma.lab.names <- vector( "list", K ) 
            for ( k in 1:K ) {
                mu.lab.names[[k]]       <- bquote( mu[.( k )] )
                sigma.lab.names[[k]]    <- bquote( sigma[.( k )] )
            }
            if (.check.grDevice() & dev ) {
                dev.new( title = paste( "Densities Feature ", rr, 
                                        " Mu", sep = "" ) )
            }
            .symmetric.Dens( mu[, rr,], mu.lab.names )
            if (.check.grDevice() & dev ) {
                dev.new( title = paste( "Densities Feature ", rr, 
                                        " Sigma", sep = "" ) )
            }
            .symmetric.Dens( sigma[, rr,], sigma.lab.names )
        }        
    }
    if ( .check.grDevice() && dev ) {
        dev.new( title = "Densities Hyperparameter C" ) 
    } 
}

".dens.Studmult"  <- function( x, dev ) 
{
    K       <- x@model@K
    r       <- x@model@r
    mu      <- x@par$mu
    sigma   <- x@par$sigma
    degf    <- x@par$df
    for ( rr in 1:r ) {
        if ( K == 1 ) {
            if (.check.grDevice() & dev ) {
                dev.new( title = paste( "Densities Feature ", rr, 
                                        " Mu", sep = "" ) )
            }
            .symmetric.Dens( mu[, rr,], list( bquote( mu ) ) )
            if (.check.grDevice() & dev ) {
                dev.new( title = paste( "Densities Feature ", rr, 
                                        " Sigma", sep = "" ) )
            }
            .symmetric.Dens( sigma[, rr,], list( bquote( sigma ) ) )
         
        } else {
            mu.lab.names    <- vector( "list", K )
            sigma.lab.names <- vector( "list", K ) 
            for ( k in 1:K ) {
                mu.lab.names[[k]]       <- bquote( mu[.( k )] )
                sigma.lab.names[[k]]    <- bquote( sigma[.( k )] )
            }
            if (.check.grDevice() & dev ) {
                dev.new( title = paste( "Densities Feature ", rr, 
                                        " Mu", sep = "" ) )
            }
            .symmetric.Dens( mu[, rr,], mu.lab.names )
            if (.check.grDevice() & dev ) {
                dev.new( title = paste( "Densities Feature ", rr, 
                                        " Sigma", sep = "" ) )
            }
            .symmetric.Dens( sigma[, rr,], sigma.lab.names )
        }
    }
    if ( K == 1 ) { 
        if (.check.grDevice() & dev ) {
            dev.new( title = "Density Degrees of Freedom" )                                      
        }
        .symmetric.Dens( degf[, rr,], list( bquote( nu ) ) ) 
    } else {
        degf.lab.names  <- vector( "list", K )
        for ( k in 1:K ) {
            degf.lab.names[[k]] <- bquote( nu[.( k )] )
        }
        if (.check.grDevice() & dev ) {
            dev.new( title = "Densities Degrees of Freedom" )                                    
        }
        .symmetric.Dens( degf[, rr,], degf.lab.names ) 
    }
}

### Plot Point Processes
### Plot Point Process Poisson: Plots the point process
### for the MCMC draws for lambda. The values are plotted
### against a random normal sample. 
".pointproc.Poisson" <- function(x, dev)
{
    K   <- x@model@K
    M   <- x@M
    if (.check.grDevice() && dev) {
        dev.new(title = "Point Process Representation (MCMC)")
    }
    y.grid  <- replicate(K, rnorm(M))
    if (median(x@par$lambda) < 1) {
        lambda  <- log(x@par$lambda)
    } else {
        lambda  <- x@par$lambda
    }
    col.grid <- gray.colors(K, start = 0, 
                            end = 0.5)
    legend.names    <- vector("list", K)
    for (k in seq(1, K)) {
        legend.names[[k]]   <- bquote(lambda[.(k)])
    }
    plot(lambda, y.grid, pch = 20, col = col.grid,
         cex = .7, cex.axis = .7, cex.lab = .7,
         main = "", ylab = "", xlab = "")
    mtext(side = 1, bquote(lambda), cex = .7, 
          cex.lab = .7, line = 3)
    legend("topright", legend = do.call(expression, 
                                        legend.names),
           col = col.grid, fill = col.grid)
}

".pointproc.Binomial" <- function(x, dev)
{
    K   <- x@model@K
    M   <- x@M
    if (.check.grDevice() && dev) {
        dev.new(title = "Point Process Representation (MCMC)")
    }
    y.grid  <- replicate(K, rnorm(M))
    p <- x@par$p
    col.grid <- gray.colors(K, start = 0, 
                            end = 0.5)
    legend.names    <- vector("list", K)
    for (k in seq(1, K)) {
        legend.names[[k]]   <- bquote(p[.(k)])
    }
    plot(p, y.grid, pch = 20, col = col.grid,
         cex = .7, cex.axis = .7, cex.lab = .7,
         main = "", ylab = "", xlab = "")
    mtext(side = 1, bquote(p), cex = .7, 
          cex.lab = .7, line = 3)
    legend("topright", legend = do.call(expression, 
                                        legend.names),
           col = col.grid, fill = col.grid)
}

### Plot sampling representation
### Plot sampling representation Poisson: Plots the sampling
### representation for Poisson parameters. Each parameter sample
### is combined with the other samples. 
".samprep.Poisson" <- function(x, dev)
{
    K       <- x@model@K
    if (K == 1) {
        warning(paste("Sampling representation is only ",
                      "available for mixture models with ",
                      "K > 1.", sep = ""))
        return(FALSE)
    }
    M       <- x@M
    n       <- min(2000, x@M)
    n.perm  <- choose(K, 2) * factorial(2)
    lambda  <- x@par$lambda
    if (.check.grDevice() && dev) {
        dev.new(title = "Sampling Representation (MCMC)")
    }
    comb    <- as.matrix(expand.grid(seq(1, K), seq(1, K)))
    comb    <- comb[which(comb[, 1] != comb[, 2]), ]
    lambda  <- lambda[seq(1, n), ]
    lambda  <- matrix(lambda[,comb], nrow = n * n.perm, ncol = 2)
    plot(lambda, col = "gray47", cex.lab = .7, cex.axis = .7,
         cex = .7, pch = 20, main = "", xlab = "", ylab = "")
    abline(0, 1, lty = 1)
    mtext(side = 1, bquote(lambda), cex = .7, cex.lab = .7,
          line = 3)
    mtext(side = 2, bquote(lambda), cex = .7, cex.lab = .7,
          line = 3)

}

".samprep.Binomial" <- function(x, dev)
{
    K       <- x@model@K
    if (K == 1) {
        warning(paste("Sampling representation is only ",
                      "available for mixture models with ",
                      "K > 1.", sep = ""))
        return(FALSE)
    }
    M       <- x@M
    n       <- min(2000, x@M)
    n.perm  <- choose(K, 2) * factorial(2)
    p       <- x@par$p
    if (.check.grDevice() && dev) {
        dev.new(title = "Sampling Representation")
    }
    comb    <- as.matrix(expand.grid(seq(1, K), seq(1, K)))
    comb    <- comb[which(comb[, 1] != comb[, 2]), ]
    p       <- p[seq(1, n), ]
    p       <- matrix(p[,comb], nrow = n * n.perm, ncol = 2)
    plot(p, col = "gray47", cex.lab = .7, cex.axis = .7,
         cex = .7, pch = 20, main = "", xlab = "", ylab = "")
    abline(0, 1, lty = 1)
    mtext(side = 1, bquote(p), cex = .7, cex.lab = .7,
          line = 3)
    mtext(side = 2, bquote(p), cex = .7, cex.lab = .7,
          line = 3)

}

### Posterior Density
### Posterior Density Poisson: Plots a contour plot of the 
### posterior density of the sampled parameters for K = 2.
".postdens.Poisson" <- function(x, dev)
{
    K   <- x@model@K
    if (K != 2) {
        warning(paste("A plot of the posterior density is ",
                      "available only for K = 2.", sep = ""))
    } else {
        M   <- x@M
        n   <- min(2000, M)
        lambda  <- x@par$lambda
        lambda  <- lambda[seq(1, n), ]
        dens    <- bkde2D(lambda, bandwidth = c(sd(lambda[, 1]),
                                                sd(lambda[, 2])))
        if (.check.grDevice() && dev) {
            dev.new(title = "Posterior Density Contour Plot (MCMC)")
        } 
        contour(dens$x1, dens$x2, dens$fhat, cex = .7, 
                cex.lab = .7, cex.axis = .7, col = "gray47", 
                main = "", xlab = "", ylab = "")
        mtext(side = 1, bquote(lambda[1]), cex = .7, 
              cex.lab = .7, line = 3)
        mtext(side = 2, bquote(lambda[2]), cex = .7,
              cex.lab = .7, line = 3)
        if (.check.grDevice() && dev) {
            dev.new(title = "Posterior Density Perspective Plot (MCMC)")
        }
        persp(dens$x1, dens$x2, dens$fhat, col = "gray65", 
              border = "gray47", theta = 55, phi = 30, 
              expand = .5, lphi = 180, ltheta = 90, 
              r = 40, d = .1, ticktype = "detailed", zlab = 
              "Density", xlab = "k = 1" , ylab = "k = 2")
    }
}

".postdens.Binomial" <- function(x, dev)
{
    K   <- x@model@K
    if (K != 2) {
        warning(paste("A plot of the posterior density is ",
                      "available only for K = 2.", sep = ""))
    } else {
        M       <- x@M
        n       <- min(2000, M)
        p       <- x@par$p
        p       <- p[seq(1, n), ]
        dens    <- bkde2D(p, bandwidth = c(sd(p[, 1]),
                                           sd(p[, 2])))
        if (.check.grDevice() && dev) {
            dev.new(title = "Posterior Density Contour Plot (MCMC)")
        } 
        contour(dens$x1, dens$x2, dens$fhat, cex = .7, 
                cex.lab = .7, cex.axis = .7, col = "gray47", 
                main = "", xlab = "", ylab = "")
        mtext(side = 1, bquote(p[1]), cex = .7, 
              cex.lab = .7, line = 3)
        mtext(side = 2, bquote(p[2]), cex = .7,
              cex.lab = .7, line = 3)
        if (.check.grDevice() && dev) {
            dev.new(title = "Posterior Density Perspective Plot (MCMC)")
        }
        persp(dens$x1, dens$x2, dens$fhat, col = "gray65", 
              border = "gray47", theta = 55, phi = 30, 
              expand = .5, lphi = 180, ltheta = 90, 
              r = 40, d = .1, ticktype = "detailed", zlab = 
              "Density", xlab = "k = 1" , ylab = "k = 2")
    }
}

### Logic
### Logic subseq: This function is used for each 
### distribution type in 'model'. It crreates a subsequence
### for the log-likelihoods. 
".subseq.Log.Fix" <- function(obj, index)
{
    obj@log$mixlik     <- matrix(obj@log$mixlik[index],
                                 nrow = obj@M, ncol = 1)
    obj@log$mixprior   <- matrix(obj@log$mixprior[index],
                                 nrow = obj@M, ncol = 1)
    return(obj)
}

### Logic subseq Poisson: This function creates a subsequence 
### MCMC Poisson parameter samples. 
".subseq.Poisson" <- function(obj, index) 
{
    if (obj@model@K == 1) {
        obj@par$lambda <- matrix(obj@par$lambda[index], 
                                    nrow = obj@M, ncol = 1)
    } else {
        obj@par$lambda <- obj@par$lambda[index,]
    }
    return(obj)
}

### 

".subseq.Binomial" <- function(obj, index) 
{
    if (obj@model@K == 1) {
        obj@par$p <- matrix(obj@par$p[index], nrow = obj@M,
                            ncol = 1)
    } else {
        obj@par$p <- obj@par$p[index,]
    }
    return(obj)
}

### 

".subseq.Normal" <- function( obj, index )
{
    if ( obj@model@K == 1 ) {
        obj@par$mu <- matrix( obj@par$mu[index], nrow = obj@M,
                             ncol = 1 )
        obj@par$sigma <- matrix( obj@par$mu[index], nrow = obj@M,
                                 ncol = 1)
    } else {
        obj@par$mu  <- obj@par$mu[index, ]
        obj@par$sigma <- obj@par$sigma[index, ]
    }
    return( obj )
}

###

".subseq.Student" <- function( obj, index ) 
{
    if ( obj@model@K == 1 ) {
        obj@par$mu      <- matrix( obj@par$mu[index], nrow = obj@M,
                                   ncol = 1 )
        obj@par$sigma   <- matrix( obj@par$sigma[index], nrow = obj@M,
                                   ncol = 1 )
        obj@par$df      <- matrix( obj@par$df[index], nrow = obj@M,
                                   ncol = 1 )
    } else {
        obj@par$mu      <- obj@par$mu[index, ]
        obj@par$sigma   <- obj@par$sigma[index, ]
        obj@par$df      <- obj@par$df[index, ]
    }
    return( obj )
}

".subseq.Normult" <- function( obj, index ) 
{
    if ( obj@model@K == 1 ) {
        obj@par$mu          <- matrix( obj@par$mu[index,], nrow = obj@M,
                                       ncol = 1 )
        obj@par$sigma       <- matrix( obj@par$sigma[index,], nrow = obj@M,
                                       ncol = 1 )        
        obj@par$sigmainv    <- matrix( obj@par$sigmainv[index, ], nrow = obj@M,
                                       ncol = 1 )
    } else {
        obj@par$mu          <- obj@par$mu[index,,]
        obj@par$sigma       <- obj@par$sigma[index,,]
        obj@par$sigmainv    <- obj@par$sigmainv[index,,]
    }
    return( obj )
}

".subseq.Studmult" <- function( obj, index ) 
{
     if ( obj@model@K == 1 ) {
        obj@par$mu          <- obj@par$mu[index,]
        obj@par$sigma       <- obj@par$sigma[index,]
        obj@par$sigmainv    <- obj@par$sigmainv[index,]
        obj@par$df          <- obj@par$df[index]
    } else {
        obj@par$mu          <- obj@par$mu[index,,]
        obj@par$sigma       <- obj@par$sigma[index,,]
        obj@par$sigmainv    <- obj@par$sigmainv[index,,]
        obj@par$df          <- obj@par$df[index,] 
    }
    return( obj )
}
### Log swapElements
### Logic swapElements Poisson: This function permutes
### the elements in the MCMC sample for Poisson 
### parameters by calling the C++-function 'swap_cc()'.
".swapElements.Poisson" <- function( obj, index )
{
    ## Rcpp::export 'swap_cc'
    obj@par$lambda <- swap_cc( obj@par$lambda, index )
    return(obj)
}

###

".swapElements.Binomial" <- function( obj, index ) 
{
    ## Rcpp::export 'swap_cc'
    obj@par$p <- swap_cc( obj@par$p, index )
    return(obj)
}

".swapElements.Exponential" <- function( obj, index )
{
    ## Rcpp::export 'swap_cc'
    obj@par$lambda  <- swap_cc( obj@par$lambda, index )
    return( obj )
}

".swapElements.Normal" <- function( obj, index )
{
    ## Rcpp::export 'swap_cc'
    obj@par$mu      <- swap_cc( obj@par$mu, index )
    obj@par$sigma   <- swap_cc( obj@par$sigma, index )
    return( obj )
}

".swapElements.Student" <- function( obj, index ) 
{    
    ## Rcpp::export 'swap_cc'
    obj@par$mu      <- swap_cc( obj@par$mu, index )
    obj@par$sigma   <- swap_cc( obj@par$sigma, index )
    obj@par$df      <- swap_cc( obj@par$df, index )
    return( obj )
}

".swapElements.Normult" <- function( obj, index) 
{
    ## Rcpp::export 'swap_3d_cc' 
    obj@par$mu          <- swap_3d_cc( obj@par$mu, index )
    obj@par$sigma       <- swap_3d_cc( obj@par$sigma, index )
    obj@par$sigmainv    <- swap_3d_cc( obj@par$sigma, index )
    return( obj )
}

".swapElements.Studmult" <- function( obj, index) 
{
    ## Rcpp::export 'swap_3d_cc' 
    obj@par$mu          <- swap_3d_cc( obj@par$mu, index )
    obj@par$sigma       <- swap_3d_cc( obj@par$sigma, index )
    obj@par$sigmainv    <- swap_3d_cc( obj@par$sigma, index ) 
    obj@par$df          <- swap_cc( obj@par$df, index )
    return( obj )
}

### Validity
### Validity subseq: The index given to 'subseq()' must 
### have dimension M x 1 and must contain logical values.
".subseq.valid.Arg" <- function(obj, index) 
{
    if (dim(index)[1] != obj@M) {
        stop("Argument 'index' has wrong dimension.")
    }
    if (typeof(index) != "logical") {
        stop("Argument 'index' must be of type 'logical'.")
    }
}

### Validity swapElements: The index given to 'swapElements()'
### must have dimension M x K. It must be of type 'integer'
### and must be in the range 1, ..., K.
".swapElements.valid.Arg" <- function(obj, index)
{
    if (dim(index)[1] != obj@M || dim(index)[2] != obj@model@K) {
        stop("Argument 'index' has wrong dimension.")
    }
    if (typeof(index) != "integer") {
        stop("Argument 'index' must be of type 'integer'.")
    }
    if (!all(index > 0) || any(index > obj@model@K)) {
        stop(paste("Elements in argument 'index' must be greater 0", 
             "and must not exceed its number of columns.", 
             sep = ""))
    }
}

### --------------------------------------------------------------
### Extract
### --------------------------------------------------------------
".extract.Normult"  <- function( obj, index )
{
    dist    <- obj@model@dist
    r       <- obj@model@r
    K       <- obj@model@K
    pars    <- sapply( obj@par, function( x ) x[index,,] )
    weight  <- as.array( obj@model@weight )
    .mcmcextract( dist = dist, K = K, r = r, par = pars, 
                  weight = weight )    
}

### --------------------------------------------------------------
### Moments
### --------------------------------------------------------------
".moments.Normult.Mcmcoutput"   <- function( obj ) 
{
    moments <- array( numeric(), dim = c( obj@M, r, ) ) 
    moments <- apply( seq( 1, obj@M ), 1, 
                      function( i ) { mm <- extract( obj, i );
                                      moms  <- moments( mm ) } )
}
