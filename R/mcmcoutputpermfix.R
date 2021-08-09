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

.mcmcoutputpermfix <- setClass("mcmcoutputpermfix",
                               contains = c("mcmcpermfix", "mcmcoutputfix"),
                               validity = function(object) 
                               {
                                   ## else: OK
                                   TRUE
                               }
)

setMethod("initialize", "mcmcoutputpermfix", 
          function(.Object, mcmcoutput, Mperm = integer(), 
                   parperm = list(), logperm = list()) 
          {
              .Object@M         <- mcmcoutput@M
              .Object@burnin    <- mcmcoutput@burnin
              .Object@ranperm   <- mcmcoutput@ranperm
              .Object@par       <- mcmcoutput@par
              .Object@log       <- mcmcoutput@log
              .Object@model     <- mcmcoutput@model
              .Object@prior     <- mcmcoutput@prior
              .Object@Mperm     <- Mperm
              .Object@parperm   <- parperm
              .Object@logperm   <- logperm
              .Object
          }
)

setMethod("show", "mcmcoutputpermfix",
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
              cat("     Mperm       :", object@Mperm, "\n")
              cat("     parperm     : List of", 
                  length(object@parperm), "\n")
              cat("     logperm     : List of",
                  length(object@logperm), "\n")
              cat("     model       : Object of class", 
                  class(object@model), "\n")
              cat("     prior       : Object of class",
                  class(object@prior), "\n")
          }
)

setMethod( "plotTraces", signature( x     = "mcmcoutputpermfix", 
                                    dev   = "ANY",
                                    lik   = "ANY",
                                    col   = "ANY" ), 
	function(x, dev = TRUE, lik = 1, col = FALSE, ...) 
    {
        dist <- x@model@dist
        if ( lik %in% c( 0, 1 ) ) {            
            if( dist == "poisson" ) {
                .permtraces.Poisson( x, dev )
            } else if ( dist == "binomial" ) {
                .permtraces.Binomial( x, dev )
            } else if ( dist == "exponential" ) {
                .permtraces.Exponential( x, dev )
            } else if ( dist == "normal" ) {
                .permtraces.Normal( x, dev )
            } else if ( dist == "student") {
                .permtraces.Student( x, dev )                 
            } else if ( dist == "normult" ) {
                .permtraces.Normult( x, dev, col )
            } else if ( dist == "studmult" ) {
                .permtraces.Studmult( x, dev, col )
            }
        }
        if ( lik %in% c( 1, 2 ) ) {
            ## log ##
            .permtraces.Log( x, dev, col )
        }
	}
)

setMethod("plotHist", signature(x   = "mcmcoutputpermfix", 
                                dev = "ANY"), 
	function(x, dev = TRUE, ...) 
    {
        dist <- x@model@dist
        if(dist == "poisson") {
            .permhist.Poisson(x, dev)
        } else if (dist == "binomial") {
            .permhist.Binomial(x, dev)
        }
    }
)

setMethod("plotDens", signature(x   = "mcmcoutputpermfix",
                                dev = "ANY"),
          function(x, dev = TRUE, ...) 
          {
              dist <- x@model@dist
              if (dist == "poisson") {
                  .permdens.Poisson(x, dev)
              } else if (dist == "binomial") {
                  .permdens.Binomial(x, dev)
              }
          }
)

setMethod("plotPointProc", signature(x      = "mcmcoutputpermfix",
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

setMethod("plotSampRep", signature(x    = "mcmcoutputpermfix",
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

setMethod("plotPostDens", signature(x   = "mcmcoutputpermfix",
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

### Traces
### Traces Poisson: Plots the traces of the Poisson parameter.
".permtraces.Poisson" <- function(x, dev)
{
    K <- x@model@K
    trace.n <- K
    if (.check.grDevice() && dev) {	
        dev.new(title = "Traceplots")
    }
    par(mfrow = c(trace.n, 1), mar = c(1, 0, 0, 0), 
        oma = c(4, 5, 4,4))
    lambda <- x@parperm$lambda
    for (k in 1:K) {
        plot(lambda[, k], type = "l", axes = F, 
             col = "gray20", xlab = "", ylab = "")
        axis(2, las = 2, cex.axis = 0.7)
        mtext(side = 2, las = 2, bquote(lambda[k = .(k)]), 
              cex = 0.6, line = 3)
	}
    axis(1)
    mtext(side = 1, "Iterations", cex = 0.7, line = 3)
}	

".permtraces.Binomial" <- function(x, dev)
{
    K <- x@model@K
    trace.n <- K
    if (.check.grDevice() && dev) {	
        dev.new(title = "Traceplots")
    }
    par(mfrow = c(trace.n, 1), mar = c(1, 0, 0, 0), 
        oma = c(4, 5, 4,4))
    p <- x@parperm$p
    for (k in 1:K) {
        plot(p[, k], type = "l", axes = F, 
             col = "gray20", xlab = "", ylab = "")
        axis(2, las = 2, cex.axis = 0.7)
        mtext(side = 2, las = 2, bquote(p[k = .(k)]), 
              cex = 0.6, line = 3)
	}
    axis(1)
    mtext(side = 1, "Iterations", cex = 0.7, line = 3)
}

### --------------------------------------------------------------------
### .permtraces.Exponential
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
".permtraces.Exponential" <- function( x, dev )
{
    K       <- x@model@K
    trace.n <- K
    if ( .check.grDevice() && dev ) {	
        dev.new( title = "Traceplots" )
    }
    par( mfrow = c( trace.n, 1 ), mar = c( 1, 0, 0, 0 ), 
         oma = c( 4, 5, 4, 4 ) )
    lambda <- x@parperm$lambda
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
### .permtraces.Normal
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
".permtraces.Normal" <- function( x, dev ) 
{
    K       <- x@model@K
    trace.n <- 2 * K
    if ( .check.grDevice() && dev ) {
        dev.new( title = "Traceplots" )        
    }
    par( mfrow = c( trace.n, 1 ), mar = c( 1, 0, 0, 0 ),
         oma = c( 4, 5, 4, 4 ) )
    mu      <- x@parperm$mu
    sigma   <- x@parperm$sigma
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
### .permtraces.Student
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
".permtraces.Student" <- function( x, dev ) 
{
    K       <- x@model@K
    trace.n <- 3 * K
    if ( .check.grDevice() && dev ) {
        dev.new( title = "Traceplots" )        
    }
    par( mfrow = c( trace.n, 1 ), mar = c( 1, 0, 0, 0 ),
         oma = c( 4, 5, 4, 4 ) )
    mu      <- x@parperm$mu
    sigma   <- x@parperm$sigma
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
### .permtraces.Normult
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
".permtraces.Normult" <- function( x, dev, col ) 
{
    K       <- x@model@K
    r       <- x@model@r
    if ( .check.grDevice() && dev ) {
        dev.new( title = "Traceplots" ) 
    }
    trace.n <- r + 2
    par( mfrow = c( trace.n, 1 ), mar = c( 1, 2, 0, 0 ),
         oma = c( 4, 5, 2, 4 ) )
    mu      <- x@parperm$mu
    sigma   <- x@parperm$sigma
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
    moms    <- permmoments_cc( x )
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
### .permtraces.Studmult
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
".permtraces.Studmult" <- function( x, dev, col ) 
{
    K       <- x@model@K
    r       <- x@model@r
    if ( .check.grDevice() && dev ) {
        dev.new( title = "Traceplots" ) 
    }
    trace.n <- r + 2
    par( mfrow = c( trace.n, 1 ), mar = c( 1, 2, 0, 0 ),
         oma = c( 4, 5, 4, 4 ) )
    mu      <- x@parperm$mu
    sigma   <- x@parperm$sigma
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
    degf    <- x@parperm$df  
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
    moms    <- permmoments_cc( x )
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

### Traces log-likelihood: Plots the traces of the log-
### likelihoods.
".permtraces.Log" <- function( x, dev, col )
{
    if( .check.grDevice() && dev ) {
        dev.new( title = "Log Likelihood Traceplots" )
    }
    if ( col ) {
        cscale  <- rainbow( 3, start = 0, end = .5 )
    } else {
        cscale  <- gray.colors( 3, start = 0, end = .15 )
    }
    par( mfrow = c( 2, 1 ), mar = c( 1, 0, 0, 0 ),
         oma = c( 4, 5, 4, 4 ) )
    mixlik <- x@logperm$mixlik
    plot( mixlik, type = "l", axes = F,
          col = cscale[3], xlab = "", ylab = "" )
    axis( 2, las = 2, cex.axis = 0.7 )
    mtext( side = 2, las = 3, "mixlik", cex = 0.6,
           line = 3 )
    mixprior <- x@logperm$mixprior
    plot( mixprior, type = "l", axes = F,
         col = cscale[2], xlab = "", ylab = "" )
    axis( 2, las = 2, cex.axis = 0.7 )
    mtext( side = 2, las = 3, "mixprior", cex = 0.6,
           line = 3 )
    axis( 1 )
    mtext( side = 1, "Iterations", cex = 0.7, line = 3 )
}

### Histograms
### Histograms Poisson: Plots histograms for all Poisson
### parameters.
".permhist.Poisson" <- function(x, dev)
{
    K <- x@model@K 
    if (.check.grDevice() && dev) {
        dev.new(title = "Histograms (permuted)")
    }
    lambda <- x@parperm$lambda
    lab.names <- vector("list", K)
    for (k in 1:K) {
        lab.names[[k]] <- bquote(lambda[.(k)])
    }
    .symmetric.Hist(lambda, lab.names)
}

".permhist.Binomial" <- function(x, dev)
{
    K <- x@model@K 
    if (.check.grDevice() && dev) {
        dev.new(title = "Histograms (permuted)")
    }
    p <- x@parperm$p
    lab.names <- vector("list", K)
    for (k in 1:K) {
        lab.names[[k]] <- bquote(p[.(k)])
    }
    .symmetric.Hist(p, lab.names)
}
### Densities
### Densities Poisson: Plots densities for all Poisson
### parameters.
".permdens.Poisson" <- function(x, dev)
{
    K <- x@model@K 
    if (.check.grDevice() && dev) {
        dev.new(title = "Densities (permuted)")
    }
    lambda <- x@parperm$lambda
    lab.names <- vector("list", K)
    for (k in 1:K) {
        lab.names[[k]] <- bquote(lambda[.(k)])
    }
    .symmetric.Dens(lambda, lab.names)
}

".permdens.Binomial" <- function(x, dev)
{
    K <- x@model@K 
    if (.check.grDevice() && dev) {
        dev.new(title = "Densities (permuted)")
    }
    p           <- x@parperm$p
    lab.names   <- vector("list", K)
    for (k in 1:K) {
        lab.names[[k]] <- bquote(p[.(k)])
    }
    .symmetric.Dens(p, lab.names)
}

### Plot Point Processes
### Plot Point Process Poisson: Plots the point process
### for the MCMC draws for lambda. The values are plotted
### against a random normal sample. 
".permpointproc.Poisson" <- function(x, dev)
{
    K   <- x@model@K
    M   <- x@M
    if (.check.grDevice() && dev) {
        dev.new("Point Process Representation (MCMC permuted)")
    }
    y.grid  <- replicate(K, rnorm(M))
    if (median(x@parperm$lambda) < 1) {
        lambda  <- log(x@parperm$lambda)
    } else {
        lambda  <- x@parperm$lambda
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

".permpointproc.Binomial" <- function(x, dev)
{
    K   <- x@model@K
    M   <- x@M
    if (.check.grDevice() && dev) {
        dev.new(title = "Point Process Representation (MCMC permuted)")
    }
    y.grid  <- replicate(K, rnorm(M))
    p       <- x@par$p
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
".permsamprep.Poisson" <- function(x, dev)
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
    lambda  <- x@parperm$lambda
    if (.check.grDevice() && dev) {
        dev.new(title = "Sampling Representation")
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

".permsamprep.Binomial" <- function(x, dev)
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
    p       <- x@parperm$p
    if (.check.grDevice() && dev) {
        dev.new(title = "Sampling Representation (MCMC permuted)")
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
".permpostdens.Poisson" <- function(x, dev)
{
    K   <- x@model@K
    if (K != 2) {
        warning(paste("A plot of the posterior density is ",
                      "available only for K = 2.", sep = ""))
    } else {
        M   <- x@M
        n   <- min(2000, M)
        lambda  <- x@parperm$lambda
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
            dev.new(title = "Posterior Density Persepctive Plot (MCMC)")
        }
        persp(dens$x1, dens$x2, dens$fhat, col = "gray65", 
              border = "gray47", theta = 55, phi = 30, 
              expand = .5, lphi = 180, ltheta = 90, 
              r = 40, d = .1, ticktype = "detailed", zlab = 
              "Density", xlab = "k = 1" , ylab = "k = 2")
    }
}

".permpostdens.Binomial" <- function(x, dev)
{
    K   <- x@model@K
    if (K != 2) {
        warning(paste("A plot of the posterior density is ",
                      "available only for K = 2.", sep = ""))
    } else {
        M   <- x@M
        n   <- min(2000, M)
        p   <- x@parperm$p
        p   <- p[seq(1, n), ]
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
            dev.new(title = "Posterior Density Persepctive Plot (MCMC)")
        }
        persp(dens$x1, dens$x2, dens$fhat, col = "gray65", 
              border = "gray47", theta = 55, phi = 30, 
              expand = .5, lphi = 180, ltheta = 90, 
              r = 40, d = .1, ticktype = "detailed", zlab = 
              "Density", xlab = "k = 1" , ylab = "k = 2")
    }
}
