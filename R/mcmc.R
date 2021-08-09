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

.mcmc <- setClass( "mcmc", 
                   representation( burnin     = "integer",
                                   M          = "integer",
                                   startpar   = "logical",
                                   storeS     = "integer",
                                   storepost  = "logical",                                 
                                   ranperm    = "logical",
                                   storeinv   = "logical" ),
                   validity = function( object ) 
                   {     
                       .valid.MCMC(object)                        
                       ## else: OK
                       TRUE
                   },
                   prototype( burnin      = integer(),
                              M           = integer(),
                              startpar    = logical(),
                              storeS      = integer(),
                              storepost   = logical(),
                              ranperm     = logical(),
                              storeinv    = logical()
                             )
)

"mcmc" <- function( burnin = 0, M = 5000,
                    startpar = TRUE, storeS = 1000, 
                    storepost = TRUE, ranperm = TRUE,
                    storeinv = TRUE ) 
{
    .mcmc( burnin = as.integer( burnin ),  
           M = as.integer( M ), startpar = startpar, 
           storeS = as.integer( storeS ), storepost = storepost, 
           ranperm = ranperm, storeinv = storeinv )
}

setMethod( "show", "mcmc", 
           function( object ) {
               cat( "Object 'mcmc'\n" )
               cat( "     class       :", class(object), "\n" )
               cat( "     burnin      :", object@burnin, "\n" )
               cat( "     M           :", object@M, "\n" )
               cat( "     startpar    :", object@startpar, "\n" )
               cat( "     storeS      :", object@storeS, "\n" )
               cat( "     storepost   :", object@storepost, "\n" )
               cat( "     ranperm     :", object@ranperm, "\n" )	
               cat( "     storeinv    :", object@storeinv, "\n" )
           }
)

## Getters ##
setMethod( "getBurnin", "mcmc", 
           function( object ) 
           {
               return( object@burnin )
           }
)

setMethod( "getM", "mcmc", 
           function( object ) 
           {
               return( object@M )
           }
)

setMethod( "getStartpar", "mcmc", 
           function( object ) 
           {
               return( object@startpar )
           }
)

setMethod( "getStoreS", "mcmc", 
          function( object ) 
          {
              return( object@storeS )
          }
)

setMethod( "getStorepost", "mcmc", 
           function( object ) 
           {         
               return( object@storepost )
           }
)

setMethod( "getRanperm", "mcmc", 
           function( object ) 
           {
               return( object@ranperm )
           }
)

## Setters ##
setReplaceMethod( "setBurnin", "mcmc", 
                  function( object, value ) 
                  {
                      object@burnin <- as.integer( value )
                      validObject( object )
                      return( object )
                  }
)

setReplaceMethod( "setM", "mcmc", 
                  function( object, value ) 
                  {
                      object@M <- as.integer( value )
                      validObject( object )
                      return( object )
                  }
)

setReplaceMethod( "setStartpar", "mcmc", 
                  function( object, value ) 
                  {
                      object@startpar <- value
                      validObject( object )
                      return( object )
                  }
)

setReplaceMethod( "setStoreS", "mcmc", 
                  function( object, value ) 
                  {
                      object@storeS <- as.integer( value )
                      validObject( object )
                      return( object )
                  }
)

setReplaceMethod( "setStorepost", "mcmc", 
                  function( object, value ) 
                  {
                      object@storepost <- value
                      validObject( object )
                      return( object )
                  }
)

setReplaceMethod( "setRanperm", "mcmc", 
                  function( object, value ) 
                  {
                      object@ranperm <- value
                      validObject( object )
                      return( object )
                  }
)

### Private functions
### These functions are not exported

### Valid mcmc: The number of burnins @burnin and the number of 
### last indicator vectors to store @storeS must be non-negative
### 'integers'. The number of MCMC draws @M must be a positive 'integer'.
".valid.MCMC" <- function( object )
{
    if ( object@burnin < as.integer( 0 ) ) {
        stop( paste( "Number of Burn-In draws in slot 'burnin' must be ",
              "nonnegative.", sep = "" ) )
    } else if ( object@M <= as.integer( 0 ) ) {
        stop( "Number of MCMC draws in slot 'M' must be positive. " )
    } else if ( object@storeS < as.integer( 0 ) ) {
        stop( paste( "Number of indicators to store in slot 'storeS' must be ",
              "nonnegative.", sep = "" ) )
    }
    if ( object@storeS > object@M ) {
        stop( paste( "Number of indicators to store in slot 'storeS' must be ",
                     "smaller or equal to the number of MCMC draws in slot 'M'.",
                     sep = "" ) )
    }
}

