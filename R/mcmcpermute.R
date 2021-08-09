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

"mcmcpermute" <- function( mcmcout, fdata = NULL, method = "kmeans" ) {
    ## Check arguments ##
    .check.arg.Mcmcpermute( mcmcout )
    match.arg( method, c( "kmeans", "Stephens1997a", "Stephens1997b" ) )
    mcmcout <- .coerce.Mcmcpermute( mcmcout )
    if ( method == "kmeans" ) {
        .kmeans.Mcmcpermute(mcmcout)
    } else if ( method == "Stephens1997a" ) {
        .stephens1997a.Mcmcpermute( mcmcout )
    } else {        
        .stephens1997b.Mcmcpermute( mcmcout, fdata )
    }
}

### Private functions.
### These functions are not exported.

### Checking
### Check arguments: Checks if the 'mcmcout' object is of 
### class 'mcmcoutput' or 'mcmcoutputperm'. If not an 
### error is thrown.
".check.arg.Mcmcpermute" <- function( obj )
{
    if ( !inherits( obj, c( "mcmcoutput", "mcmcoutputperm" ) ) ) {
        stop( paste( "Unkown argument: Argument 1 must inherit either from",
                     "class 'mcmcoutput' or from class 'mcmcoutputperm'.",
                     sep = "" ) )
    }
    if ( obj@model@indicfix ) {
        warning( paste( "Slot 'indicfix' of 'model' object is ",
                        "set to TRUE. For a model with fixed ",
                        "indicators no permutations can be done.",
                        sep = "" ) )
        return( obj )
    }
    if ( obj@model@K == 1 ) {
        warning( paste( "Slot 'K' of model object is set to one. ",
                        "For a model with only one component no ",
                        "permutations can be done.", sep = "" ) )
        return( obj )
    }
}

### Coercing
### Coerces any 'mcmcoutputperm' object to its corresponding
### 'mcmcoutput' object.
".coerce.Mcmcpermute" <- function(obj)
{
    ## If object is of class 'mcmcoutputperm' coerce it 
    ## to an object of class 'mcmcoutput' 
    if(inherits(obj, "mcmcoutputperm")) {
        if (class(obj) == "mcmcoutputpermfix") {
            obj <- as(obj, "mcmcoutputfix")
        } else if (class(obj) == "mcmcoutputpermfixhier") {
            obj <- as(obj, "mcmcoutputfixhier") 
        } else if (class(obj) == "mcmcoutputpermfixpost") {
            obj <- as(obj, "mcmcoutputfixpost")
        } else if (class(obj) == "mcmcoutputpermfixhierpost") {
            obj <- as(obj, "mcmcoutputfixhierpost")
        } else if (class(obj) == "mcmcoutputpermbase") {
            obj <- as(obj, "mcmcoutputbase")
        } else if (class(obj) == "mcmcoutputpermhier") {
            obj <- as(obj, "mcmcoutputhier")
        } else if (class(obj) == "mcmcoutputpermpost") {
            obj <- as(obj, "mcmcoutputpost")
        } else {
            obj <- as(obj, "mcmcoutputhierpost")
        }
    }
    return(obj)
}

### Permutation
### Kmeans: Permutes the parameters regarding a cluster
### of the vectorized parameter matrix. This can lead
### to a reduced number of MCMC draws, as values of the
### same draw could be assigned to the same component
### and are deleted from the sample.
### If no permutation is possible a warning is thrown.
### See .process.output.empty.Mcmcpermute().
".kmeans.Mcmcpermute" <- function( obj ) 
{
    K           <- obj@model@K
    M           <- obj@M
    r           <- obj@model@r
    dist        <- obj@model@dist
    ## Calculate maximum a posterior estimate (MAP)
    map.index   <- .map.Mcmcestimate( obj )
    map         <- .extract.Mcmcestimate( obj, map.index )    
    if ( dist == "poisson" ) {
        clust.par       <- sqrt( obj@par$lambda )
        clust.par       <- as.vector( clust.par )
        clust.center    <- sqrt( map$par$lambda )
    } else if ( dist == "binomial" ) {
        clust.par       <- sqrt( obj@par$p )
        clust.par       <- as.vector( clust.par )
        clust.center    <- sqrt( map$par$p )
    } else if ( dist == "normal" ) {
        clust.par       <- obj@par 
        clust.par$sigma <- sqrt( clust.par$sigma )
        clust.par       <- cbind( as.vector( clust.par$mu ), as.vector( clust.par$sigma ) )
        clust.center    <- cbind( map$par$mu, sqrt( map$par$sigma ) )
    } else if ( dist == "student" ) {
        clust.par       <- obj@par
        clust.par$sigma <- sqrt( clust.par$sigma )
        clust.par       <- cbind( as.vector( clust.par$mu ), as.vector( clust.par$sigma ),
                                  as.vector( clust.par$df ) )
        clust.center    <- cbind( map$par$mu, sqrt( map$par$sigma ), map$par$df ) 
    } else if ( dist == "normult" ) {
        clust.par       <- obj@par
        indexdiag       <- diag( qinmatr( seq( 1, r * ( r + 1 ) / 2 ) ) )
        indexdiag2      <- diag( matrix( seq( 1, r^2 ), nrow = r, byrow = TRUE ) ) 
        clust.par.sigma <- matrix( aperm( obj@par$sigma[, indexdiag, ], c( 3, 1, 2 ) ),
                                   nrow = K * M )
        clust.par.mu    <- matrix( aperm( obj@par$mu, c( 3, 1, 2 ) ),
                                   nrow = K * M )
        clust.par       <- cbind( clust.par.mu, clust.par.sigma )
        clust.center    <- cbind( t( map$par$mu ), t( qincolmult( map$par$sigma )[indexdiag,] ) ) 
    } else if ( dist == "studmult" ) {
        clust.par       <- obj@par
        indexdiag       <- diag( qinmatr( seq( 1, r * ( r + 1 ) / 2 ) ) )
        indexdiag2      <- diag( matrix( seq( 1, r^2 ), nrow = r, byrow = TRUE ) ) 
        clust.par.sigma <- matrix( aperm( obj@par$sigma[, indexdiag, ], c( 3, 1, 2 ) ),
                                   nrow = K * M )
        clust.par.mu    <- matrix( aperm( obj@par$mu, c( 3, 1, 2 ) ),
                                   nrow = K * M )
        clust.par       <- cbind( clust.par.mu, clust.par.sigma )
        clust.center    <- cbind( t( map$par$mu ), t( qincolmult( map$par$sigma )[indexdiag,] ) )   
    }
    ## Apply unsupervised k-means clustering to parameters
    if ( dist %in% c( "poisson", "binomial", "exponential" ) ) {
        result.clust    <- kmeans( clust.par, centers = as.vector( clust.center ) )
    } else {
        result.clust    <- kmeans( clust.par, centers = clust.center )
    }
    ## Parameters have been stacked vertically into a vector. Reorder.
    perm.index      <- array( result.clust$clust, dim = c( M, K ) )
    ## Check if each cluster has been hit. 
    comp.index      <- as.array( matrix( seq( 1:K ), nrow = M, ncol = K, 
                                         byrow = TRUE ) )
    keep.index      <- ( t( apply( perm.index, 1, sort, FALSE ) )  
                         == comp.index )
    is.perm         <- array( apply( keep.index, 1, all ) ) 
    nonperm         <- sum( !is.perm )
    if ( nonperm < M ) {
        ## Create a subsequence of the MCMC output
        obj.subseq <- subseq( obj, is.perm )
        ## Apply permutation suggested by kmeans clustering
        obj.swap <- swapElements( obj.subseq, perm.index[is.perm,] )
        ## Create 'mcmcoutputperm' objects
        .process.output.Mcmcpermute( obj, obj.swap, method = "kmeans" )
    } else {
        .process.output.empty.Mcmcpermute( obj, method = "kmeans" )
    }
}

### Stephens1997a calling function: Calls the appropriate
### algorithm to perform a relabeling following Stephens (1997a).
### The algorithm maximizes in each iteration the logsum of the 
### posterior likelihoods of all the parameter draws and chooses
### afterwards the best permutation of each parameter draw.
### If the log value does not change anymore, convergence
### is reached.
### If no permutation is possible, a warning is thrown.
### See .process.output.empty.Mcmcpermute().
".stephens1997a.Mcmcpermute" <- function(obj) 
{
    dist <- obj@model@dist
    ## Apply Stephens1997a relabeling algorithm
    if (dist == "poisson" || dist == "exponential" ) {
        index <- .stephens1997a.poisson.Mcmcpermute(obj)
    }
    if (dist == "binomial") {
        index <- .stephens1997a.binomial.Mcmcpermute(obj)
    } 
    ## Create 'mcmcoutputperm' objects
    startidx    <- matrix(seq(1, obj@model@K), nrow = obj@M, 
                          ncol = obj@model@K, byrow = TRUE)
    if (!identical(startidx, index)) {
        obj.swap    <- swapElements(obj, index)
        .process.output.Mcmcpermute(obj, obj.swap, method = "Stephens1997a")
    } else {
        .process.output.empty.Mcmcpermute(obj, method = "Stephens1997a")
    }                   
}

### Stephens1997a Poisson: Specific algorithm for Stephens 
### relabeling of Poisson mixtures. In this case a bounded
### Nelder-Mead algorithm is used from the 'dfoptim' package.
".stephens1997a.poisson.Mcmcpermute" <- function(obj)
{
    M                   <- obj@M
    K                   <- obj@model@K
    w.mean              <- apply(obj@post$weight, 2, mean)
    a.mean              <- apply(obj@post$par$a, 2, mean)
    b.mean              <- apply(obj@post$par$b, 2, mean)
    startpar            <- c(w.mean, a.mean, b.mean)   
    lambda              <- obj@par$lambda
    weight              <- obj@weight
    index               <- array(integer(), dim = c(M, K))
    storage.mode(index) <- "integer"
    index.out           <- matrix(seq(1, K), nrow = M, ncol = K, byrow = TRUE)
    storage.mode(index.out) <- "integer"
    perm                <- as.matrix(expand.grid(seq(1, K), seq(1, K)))
    ind                 <- apply(perm, 1, function(x) all(x == x[1]))
    perm                <- perm[!ind, ]
    storage.mode(perm)  <- "integer"
    stephens1997a_poisson_cc(lambda, weight, startpar, perm)
}

".stephens1997a.binomial.Mcmcpermute" <- function(obj)
{
    M                   <- obj@M
    K                   <- obj@model@K
    w.mean              <- apply(obj@post$weight, 2, mean)
    a.mean              <- apply(obj@post$par$a, 2, mean)
    b.mean              <- apply(obj@post$par$b, 2, mean)
    startpar            <- c(w.mean, a.mean, b.mean)   
    p                   <- obj@par$p
    weight              <- obj@weight
    index               <- array(integer(), dim = c(M, K))
    storage.mode(index) <- "integer"
    index.out           <- matrix(seq(1, K), nrow = M, ncol = K, byrow = TRUE)
    storage.mode(index.out) <- "integer"
    perm                <- as.matrix(expand.grid(seq(1, K), seq(1, K)))
    ind                 <- apply(perm, 1, function(x) all(x == x[1]))
    perm                <- perm[!ind, ]
    storage.mode(perm)  <- "integer"
    stephens1997a_poisson_cc(p, weight, startpar, perm)
}

### Stephens1997b calling function: Calls the approrpiate
### algorithm to perform a relabeling following Stephens (1997b) to 
### the sample draws and the data. The algorithm computes
### the classification probability matrices for each parameter
### draw from the MCMC sample. It computes then the best estimate
### of the classification probability matrix and the corresponding
### Kullback-Leibler distance of each matrix to this estimate. 
### For each parameter an each component of the estimate the
### 'cost' matrix is then minimized and the assignment indicates
### the assignment of a parameter draw to its label.
### This function needs an 'fdata' object and checks within
### if it is valid in regard to the 'model' object carried
### by the 'mcmcoutput' object. 
### If no permutation is possible, a warning is thrown.
".stephens1997b.Mcmcpermute" <- function( obj, fdata.obj ) 
{
    .check.fdata.model.Mcmcstart( fdata.obj, obj@model )
    dist <- obj@model@dist
    if ( dist == "poisson" ) {
        index <- .stephens1997b.poisson.Mcmcpermute( obj, fdata.obj )
    } else if ( dist == "binomial" ) {
        index <- .stephens1997b.binomial.Mcmcpermute( obj, fdata.obj ) 
    } else if ( dist == "exponential" ) {
        index <- .stephens1997b.exponential.Mcmcmpermute( obj, fdata.obj )
    }
    ## Create 'mcmcoutputperm' objects. 
    startidx    <- matrix( seq( 1, obj@model@K ), nrow = obj@M, 
                           ncol = obj@model@K, byrow = TRUE )
    if ( any( startidx != index ) ) {
        obj.swap    <- swapElements( obj, index )
        .process.output.Mcmcpermute( obj, obj.swap, method = "Stephens1997b" )
    } else {
        .process.output.empty.Mcmcpermute( obj, method = "Stephens1997b" )
    }
}

".stephens1997b.poisson.Mcmcpermute" <- function( obj, fdata.obj )  
{
    stephens1997b_poisson_cc( fdata.obj@y, obj@par$lambda,
                                       obj@weight )
}

".stephens1997b.binomial.Mcmcpermute" <- function( obj, fdata.obj )
{
    stephens1997b_binomial_cc( fdata.obj@y, fdata.obj@T, obj@par$p,
                               obj@weight )
}

".stephens1997b.exponential.Mcmcpermute" <- function( obj, fdata.obj ) 
{
    stephens1997b_exponential_cc( fdata.obj@y, obj@par$lambda, 
                                  obj@ weight )
}

".process.output.Mcmcpermute" <- function( obj, obj.swap, method )
{
    ## Create 'mcmcoutputperm' objects ##
    if ( class( obj ) == "mcmcoutputfix" ) {
        .mcmcoutputpermfix( obj, 
                            Mperm        = obj.swap@M,
                            parperm      = obj.swap@par,
                            logperm      = obj.swap@log )
    } else if ( class( obj ) == "mcmcoutputfixhier" ) {
        .mcmcoutputpermfixhier( obj, 
                                Mperm        = obj.swap@M,
                                parperm      = obj.swap@par,
                                logperm      = obj.swap@log )
    } else if (class(obj) == "mcmcoutputfixpost") {
        .mcmcoutputpermfixpost(obj,
                               Mperm        = obj.swap@M,
                               parperm      = obj.swap@par,
                               logperm      = obj.swap@log,
                               postperm     = obj.swap@post)
    } else if (class(obj) == "mcmcoutputfixhierpost") {
        .mcmcoutputpermfixhierpost(obj,
                                   Mperm        = obj.swap@M,
                                   parperm      = obj.swap@par,
                                   logperm      = obj.swap@log,
                                   postperm     = obj.swap@post)
    } else if (class(obj) == "mcmcoutputbase") {
        .mcmcoutputpermbase(obj,
                            Mperm        = obj.swap@M,
                            parperm      = obj.swap@par,       
                            relabel      = method,
                            weightperm   = obj.swap@weight,
                            logperm      = obj.swap@log,
                            entropyperm  = obj.swap@entropy,
                            STperm       = obj.swap@ST,
                            Sperm        = obj.swap@S,
                            NKperm       = obj.swap@NK)          
    } else if (class(obj) == "mcmcoutputhier") {
        .mcmcoutputpermhier(obj,
                            Mperm        = obj.swap@M,
                            parperm      = obj.swap@par,
                            relabel      = method,
                            weightperm   = obj.swap@weight,
                            logperm      = obj.swap@log,
                            entropyperm  = obj.swap@entropy,
                            STperm       = obj.swap@ST,
                            Sperm        = obj.swap@S,
                            NKperm       = obj.swap@NK)
    } else if (class(obj) == "mcmcoutputpost") {
        .mcmcoutputpermpost(obj,
                            Mperm        = obj.swap@M,
                            parperm      = obj.swap@par,
                            relabel      = method,
                            weightperm   = obj.swap@weight,
                            logperm      = obj.swap@log,
                            postperm     = obj.swap@post,
                            entropyperm  = obj.swap@entropy,
                            STperm       = obj.swap@ST,
                            Sperm        = obj.swap@S,
                            NKperm       = obj.swap@NK)           
    } else {
        .mcmcoutputpermhierpost(obj,
                                Mperm        = obj.swap@M,
                                parperm      = obj.swap@par,
                                relabel      = method,
                                weightperm   = obj.swap@weight,
                                logperm      = obj.swap@log,
                                postperm     = obj.swap@post,
                                entropyperm  = obj.swap@entropy,
                                STperm       = obj.swap@ST,
                                Sperm        = obj.swap@S,
                                NKperm       = obj.swap@NK)           
    }
}

".process.output.empty.Mcmcpermute" <- function(obj, method)
{
    warning(paste("Not a single draw is a permutation in the ",
                  "function 'mcmcpermute()'.", sep = ""))
    ## Create 'mcmcoutputperm' objects ##
    if (class(obj) == "mcmcoutputfix") {
        .mcmcoutputpermfix(obj, Mperm = as.integer(0)) 
    } else if (class(obj) == "mcmcoutputfixhier") {
        .mcmcoutputpermfixhier(obj, Mperm = as.integer(0)) 
    } else if (class(obj) == "mcmcoutputfixpost") {
        .mcmcoutputpermfixpost(obj, Mperm = as.integer(0))
    } else if (class(obj) == "mcmcoutputfixhierpost") {
        .mcmcoutputpermfixhierpost(obj, Mperm = as.integer(0))
    } else if (class(obj) == "mcmcoutputbase") {
        .mcmcoutputpermbase(obj, Mperm = as.integer(0), relabel = method)
    } else if (class(obj) == "mcmcoutputhier") {
        .mcmcoutputpermhier(obj, Mperm = as.integer(0), relabel = method)
    } else if (class(obj) == "mcmcoutputpost") {
        .mcmcoutputpermpost(obj, Mperm = as.integer(0), relabel = method) 
    } else {
        .mcmcoutputpermhierpost(obj, Mperm = as.integer(0), relabel = method)
    }
}


