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

"mcmcstart" <- function( fdata, model, varargin ) 
{
    ## Check arguments 
    .check.fdata.model.Mcmcstart( fdata, model )
	## Check if mcmc object was given in arguments
	if ( nargs() == 2 ) {
		mcmc <- mcmc()
	}		
	else {
        .check.mcmc.Mcmcstart( varargin )
		mcmc <- varargin
	}
    K       <- model@K
    dist    <- model@dist
    ## If @startpar is
    ## TRUE (default):      Start by sampling the parameters
    ##                          -> it needs starting indicators S
    ##                      In this case staring indicators are generated 
    ##                      by kmeans-clustering independently of
    ##                      the model.
    ## FALSE:               Start by sampling the indicators
    ##                          -> it needs starting parameters model@par
    ##                      In this case starting parameters are generated
    ##                      dependent on the data in @y of 'fdata' and the 
    ##                      the model.
    ## If the model has fixed indicators (@indicfix = TRUE), no indicators
    ## are generated.
    if ( !model@indicfix ) {
        if ( mcmc@startpar ) {
            if ( K > 1 ) {
                fdata <- .indicators.Mcmcstart( fdata, model )
                if ( model@dist %in% c( "student", "studmult" ) ) {
                    # As an independent prior samples conditional on the means
                    # starting means are also needed. This can be simplified 
                    # in the future by defining only starting values for the mean 
                    # and not the variances.
                    model   <- .parameters.Mcmcstart( fdata, model, mcmc )
                    model@par$sigma <- NULL
                }                
            }
        } else {
            model <- .parameters.Mcmcstart( fdata, model, mcmc )
        }
    } else {
        warning( paste( "Slot 'indicfix' of 'model' object is ",
                        "set to TRUE. 'mcmcstart()' does not ",
                        "generate indicators nor starting parameters ",
                        "for models with fixed indicators.", sep = "" ) )
    }
	obj.list <- list( fdata = fdata, model = model, mcmc = mcmc )
	return( obj.list )
}

### Private functions.
### These functions are not exported.

### Checking
### Check fdata/model: 'fdata' must be a valid 'fdata' object and 'model'
### must be a valid 'model' object. Furthermore, 'fdata' must have a 
### non-empty data slot @y. 
### If the distributions in 'model' do not correspond to the dimensions
### @r in 'fdata' an error is thrown.
".check.fdata.model.Mcmcstart" <- function( fdata.obj, model.obj ) 
{
    .valid.Fdata( fdata.obj )    
    .valid.dist.Model( model.obj )
    .valid.K.Model( model.obj )
    .valid.r.Model( model.obj )
    .valid.T.Model( model.obj )
    hasY( fdata.obj, verbose = TRUE )
    if ( fdata.obj@r > 1 && model.obj@dist %in% .get.univ.Model() ) {
        stop( paste( "Wrong specification of slot 'r' in 'fdata' object. ",
                     "Univariate distribution in slot 'dist' of 'model' ",
                     "object but dimension in slot 'r' of 'fdata' object ",
                     "greater 1.", sep = "" ) )
    } else if ( fdata.obj@r < 2 && model.obj@dist %in% .get.multiv.Model() ) {
        stop( paste( "Wrong specification of slot 'r' ind 'fdata' object ",
                     "Multivariate distribution in slot 'dist' if 'model' ",
                     "object but dimension in slot 'r' of 'fdata' object ",
                     "less than two.", sep = "" ) )
    }
}

### Check varargin: Argument 'varargin' must be an object of class 'mcmc'.
".check.mcmc.Mcmcstart" <- function(mcmc.obj) 
{
    if (class(mcmc.obj) != "mcmc") {
        stop(paste("Wrong argument. 'mcmc' must be an object of class ",
                   "'mcmc'.", sep = ""))
    }
    .valid.MCMC(mcmc.obj)
}
### Logic
### Logic parameters: Generates starting parameters for @dist in 
### 'model.obj'. Returns a 'model' object with starting parameters 
### in @par. 
".parameters.Mcmcstart" <- function( fdata.obj, model.obj, mcmc.obj )
{
    K       <- model.obj@K
    dist    <- model.obj@dist
    ## Check if model object for student-t distributions has 
    ## a parameter 'df'.
    if ( dist %in% c( "student", "studmult" ) ) {
        model.obj   <- .mcmcstart.Student.Df(model.obj)
    }
    ## Check if weights have been already initialized
    if ( K > 1 ) {
        if ( model.obj@indicmod == "multinomial" ) { 
             model.obj <- .parameters.multinomial.Mcmcstart( model.obj )
        } ## else: Markov model, implemented later.
    }
    if ( dist %in% c( "poisson", "cond.poisson" ) ) {
        .parameters.poisson.Mcmcstart( fdata.obj, model.obj )
    } else if ( dist == "exponential" ) {
        .parameters.exponential.Mcmcstart( fdata.obj, model.obj, mcmc.obj )
    } else if ( dist == "binomial" ) {
        .parameters.binomial.Mcmcstart(fdata.obj, model.obj)
    } else if ( dist %in% c( "normal", "student" ) ) {
        .mcmcstart.Norstud.Model(fdata.obj, model.obj, mcmc.obj)
    } else if ( dist %in% c( "normult", "studmult" ) ) {
        .mcmcstart.Norstudmult.Model( fdata.obj, model.obj, mcmc.obj )
    }
}

".mcmcstart.Exp" <- function(data.obj) 
{
    r           <- data.obj@r
    N           <- data.obj@N
    has.exp     <- (length(data.obj@exp) > 0) 
    if (has.exp) {
        if (data.obj@bycolumn) {
            if (nrow(data.obj@exp) != N && nrow(data.obj@exp) != 1) {
                stop(paste("Dimension of slot 'exp' of 'data' object",
                           "does not match dimension of slot 'y' of",
                           "'data' object."), 
                     sep = "")
            } else if (nrow(data.obj@exp) == N) {
                exp <- data.obj@exp
            } else {
                exp <- matrix(data.obj@exp[1, 1], nrow = N, ncol = 1)
            }
        } else {
            if (ncol(data.obj@exp) != N && ncol(data.obj@exp) != 1) {
                stop(paste("Dimension of slot 'exp' of 'data' object",
                           "does not match dimension of slot 'y' of",
                           "'data' object."),
                     sep = "")
            } else if (ncol(data.obj@exp) == N) {
                exp <- t(data.obj@exp)
            } else {
                exp <- matrix(data.obj@exp[1, 1], nrow = N, ncol = 1)
            }           
        }
    } else {
        exp <- matrix(1, nrow = N, ncol = 1)
    }
    return(exp)
}

".parameters.multinomial.Mcmcstart" <- function(model.obj)
{
    K   <- model.obj@K
    if (!hasWeight(model.obj)) {
        model.obj@weight <- matrix(1/K, nrow = 1, ncol = K)
    }
    return(model.obj)
}

".parameters.poisson.Mcmcstart" <- function(fdata.obj, model.obj) 
{
    K           <- model.obj@K
    datam       <- getColY(fdata.obj)
    if (!hasPar(model.obj)) {
        if (K == 1) { 
            pm <- max(mean(datam/exp, na.rm = TRUE), 0.1)
            pm <- array(pm, dim = c(1, K))
        } else { ## K > 1
            if (hasExp(fdata.obj)) {
                expos   <- getColExp(fdata.obj)
                pm      <- (mean(datam/expos, na.rm = TRUE)) * exp(runif(K))
            } else {
                pm      <- (mean(datam, na.rm = TRUE)) * exp(runif(K))
            }
            pm <- pmax(pm, 0.1)
        }
        model.obj@par <- list(lambda = pm)
    }
    return(model.obj)
}

".parameters.exponential.Mcmcstart" <- function( fdata.obj, model.obj,
                                                 mcmc.obj )
{
    if ( !hasPar( model.obj ) ) { 
        datam   <- getColY( fdata.obj )
        K       <- model.obj@K
        if (K == 1) {
            pm <- 1/mean( datam, na.rm = TRUE )
        } else { ## K > 1
            pm <- exp( runif( K ) )/mean( datam, na.rm = TRUE )
        }
        model.obj@par <- list( lambda = pm )    
    }
    return( model.obj )
}

".parameters.binomial.Mcmcstart" <- function(fdata.obj, model.obj)                                
{
    if (!hasPar(model.obj) && hasT(fdata.obj, verbose = TRUE)) {
        datam   <- getColY(fdata.obj)
        K       <- model.obj@K           
        if (K == 1) {
            pm  <- mean(datam/fdata.obj@T, na.rm = TRUE)
            pm  <- pmin(pmax(pm, 0.1),0.9)
        } else { ## K > 1
            pm  <- mean(datam/fdata.obj@T, na.rm = TRUE) * exp(.2 * runif(K))
            pm  <- pmin(pmax(pm, 0.1), 0.9)
        }
        model.obj@par <- list(p = pm)
    }
    return(model.obj)
}

".mcmcstart.Norstud.Model" <- function( fdata.obj, model.obj,
                                        mcmc.obj )
{
    datam       <- getColY( fdata.obj )
    K           <- model.obj@K
    has.par     <- ( length( model.obj@par ) > 0 )
    start.mu    <- FALSE
    start.sigma <- FALSE
    if( !has.par ) {
        start.mu    <- TRUE
        start.sigma <- TRUE
    } else { ## has already parameters 
        start.mu      <- !"mu" %in% names( model.obj@par )
        start.sigma   <- !"sigma" %in% names( model.obj@par ) 
    }		
    if( start.mu ) {
        if( K == 1 ) {
            pm <- mean( datam, na.rm = TRUE )  
        } else { ## K > 1
            pm <- mean( datam, na.rm = TRUE ) + 
            sd( datam, na.rm = TRUE ) * runif( K )  
            pm <- matrix( pm, nrow = 1, ncol = K )
        }
        if( start.sigma ) {
            model.obj@par       <- list( mu = pm )
        } else {
            model.obj@par$mu    <- pm
        }
    }
    if( start.sigma ) {			
        pm                  <- sd( datam, na.rm = TRUE )
        pm                  <- matrix( pm, nrow = 1, ncol = K )
        model.obj@par$sigma <- pm
    }
    
    return( model.obj )
}

".mcmcstart.Norstudmult.Model" <- function( fdata.obj, model.obj,
                                            mcmc.obj )
{
    K           <- model.obj@K
    r           <- model.obj@r
    has.par     <- ( length( model.obj@par ) > 0 )
    datam       <- getColY( fdata.obj )
    ## Check if parameters are already provided ##
    start.mu    <- FALSE
    start.sigma <- FALSE
    if ( !has.par ) {
        start.mu    <- TRUE
        start.sigma <- TRUE
    } else {
        has.mu      <- "mu" %in% names( model.obj@par )
        has.sigma   <- "sigma" %in% names( model.obj@par )
        if ( !has.mu ) {
            start.mu    <- TRUE
        }
        if ( !has.sigma ) {
            start.sigma <- TRUE
        }
    }			
    cov.m <- cov( datam )	
    if (start.mu) {
        if (K == 1) {
            pm.mu   <- apply(datam, 2, mean, na.rm = TRUE)
            pm.mu   <- array(pm.mu, dim = c(1, K))
        }
        else { ## K > 1
            mean    <- apply(datam, 2, mean, na.rm = TRUE)
            pm.mu   <- matrix(0, nrow = r, ncol = K)
            for(i in 1:K) {
                pm.mu[,i] <-  matrix(mean) + t(chol(cov.m)) %*% matrix(runif(K))
            }
        }
        if (!has.par) {
            model.obj@par <- list(mu = pm.mu)
        }
        else {
            model.obj@par$mu <- pm.mu
        }
    }
    if (start.sigma) {
        model.obj@par$sigma <- array(cov.m, dim = c(r, r, K))
    }
    return(model.obj)
}

".mcmcstart.Student.Df" <- function( model.obj ) 
{
    K           <- model.obj@K
    has.par     <- ( length( model.obj@par ) > 0 )
    if ( has.par ) {
        has.df  <- "df" %in% names( model.obj@par )
        if ( !has.df ) {	
            model.obj@par$df <- array( 10, dim = c( 1, K ) )
            validObject( model.obj )	
        }			
    } else {
        model.obj@par <- list( df = array( 10, dim = c( 1, K ) ) )
    }
    return( model.obj )
}

### Logic indicators: Returns an 'fdata' object with generated 
### indicators. 
".indicators.Mcmcstart" <- function(fdata.obj, model.obj)
{
    dist    <- model.obj@dist
    if ( dist %in% c( "poisson", "cond.poisson", "exponential" ) ) {
        .indicators.poisson.Mcmcstart(fdata.obj, model.obj)
    } else if ( dist == "binomial" ) {
        .indicators.binomial.Mcmcstart(fdata.obj, model.obj)
    } else if( dist %in% c( "normal", "normult", "student", "studmult" ) ) {
        .mcmcstart.Ind.Norstud( fdata.obj, model.obj )
    } 
}

### Logic indicators for Poisson: If it is started by sampling
### the parameters a simple kmeans-clustering is performed
### to find initial indicators. If indicators are already
### in slot @S of the 'fdata' object, the 'fdata' object is
### immediately returned. 
".indicators.poisson.Mcmcstart" <- function(fdata.obj, model.obj)
{
    K           <- model.obj@K
    if ( !hasS( fdata.obj ) ) {
        datam   <- getColY( fdata.obj ) 
        S       <- matrix(kmeans( datam^.5, centers = K,
                                  nstart = K )$cluster )
        if ( fdata.obj@bycolumn ) {
            fdata.obj@S  <- S 
        } else {
            fdata.obj@S  <- t( S )
        }
    }
    return( fdata.obj )
}

".indicators.binomial.Mcmcstart" <- function( fdata.obj, model.obj ) 
{
    if ( !hasS( fdata.obj ) ) {
        K           <- model.obj@K
        datam       <- getColY( fdata.obj )
        if( ( max( datam ) - min( datam ) ) > 2 * K ) {
            ## use k-means to determine a starting classification
            if ( fdata.obj@bycolumn ) {
                fdata.obj@S <- as.matrix( kmeans( datam^.5, 
                                                  centers = K, 
                                                  nstart = K )$cluster )
            } else {
                fdata.obj@S <- t( as.matrix( kmeans( datam^.5, 
                                                     centers = K,
                                                     nstart = K )$cluster ) ) 
            }
        } else {
            ## random classification
            N           <- fdata.obj@N
            if ( fdata.obj@bycolumn ) {
                fdata.obj@S  <- as.matrix( sample( c( 1:K ), N, 
                                                   replace = TRUE) )
            } else {
                fdata.obj@S  <- t( as.matrix( sample( c( 1:K ), N, 
                                                      replace = TRUE ) ) )
            }
        }
    } 
    if ( !hasT( fdata.obj ) || length( fdata.obj@T ) != fdata.obj@N ) {
        if ( fdata.obj@bycolumn ) {
            fdata.obj@T <- matrix( as.integer( 1 ), nrow = fdata.obj@N, ncol = 1 )
        } else {
            fdata.obj@T <- matrix( as.integer( 1 ), nrow = 1, ncol = fdata.obj@N )
        }
    }
    return( fdata.obj )
}

".mcmcstart.Ind.Norstud" <- function( data.obj, model.obj ) 
{
    if ( !hasS( data.obj ) ) {
        K           <- model.obj@K
        datam       <- getColY( data.obj )
        if ( data.obj@bycolumn ) {
            data.obj@S  <- as.matrix( kmeans( datam^.5, 
                                              centers = K, 
                                              nstart = K )$cluster )
        } else {
            data.obj@S  <- t( as.matrix( kmeans( datam^.5, 
                                                 centers = K,
                                                 nstart = K )$cluster ) )
        }
    }
    return( data.obj )
}
