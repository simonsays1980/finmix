"qinmatr" <- function( q ) 
{
    if ( length( dim( q ) ) > 0 ) {
        stop( paste( "The argument 'q' has to be an object of dimension 1 x r or r x 1.",
                     sep = ""),
              call. = TRUE )        
    } 
    r   <- -.5 + sqrt( .25 + 2 * length( q ) )
    tmp <- matrix( numeric(), nrow = r, ncol = r )
    tmp[upper.tri( tmp, diag = TRUE )] <- q
    tmp[lower.tri( tmp )] <- t( tmp[upper.tri( tmp )] )
    return( tmp )
}

"qinmatrmult" <- function( m ) 
{
    r   <- -.5 + sqrt( .25 + 2 * nrow( m ) )
    tmp.array   <- array( numeric(), dim = c( r, r, ncol( m ) ) )
    for (k in 1:ncol( m ) ) {
        tmp.array[,, k]  <- qinmatr( m[, k] )
    }
    return( tmp.array )
}

"qincol" <- function( m ) 
{
    r       <- ncol( m )
    index   <- 0 
    s       <- r * ( r + 1 ) / 2
    qcol    <- vector( "numeric", s )
    for (rr in 1:r ) {
        qcol[( index + 1 ) : ( index + rr )] <- m[1:rr, rr]                  
        index   <- index + rr
    }
    return( qcol )
}

"qincolmult"    <- function( a ) 
{
    r   <- dim( a )[1] 
    K   <- dim( a )[3]
    s   <- r * ( r + 1 ) / 2
    tmp.mat     <- matrix( numeric(), nrow = s, ncol = K)
    for ( k in 1:K ) {
        tmp.mat[, k]    <- qincol( a[,, k] )
    }
    return( tmp.mat )
}

