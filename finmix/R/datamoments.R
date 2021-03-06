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

## 'datamoments' is a virtual class from which the corresponding 
## datamoments for 'continuous' and 'discrete' inherit	 		 
.datamoments <- setClass("datamoments",
                         representation(mean            = "numeric",
                                        var             = "matrix",
                                        fdata           = "fdata",
                                        "VIRTUAL"
                                        )
)

## mutual constructor for all type of datamoments ##
"datamoments" <- function(value = fdata()) 
{
        hasY(value, verbose = TRUE)
        if (value@type == "continuous") {
            tryCatch(
			    {
                    .Object <- .cdatamoments(value = value)
                }, 
                error = function( e ) 
                {
                    e$call      <- call( "datamoments" )
                    e$message   <- paste0( e$message, "\nIs your data possibly discrete valued? ",
                                           "In this case change the data type in the",
                                           " 'fdata' object to 'discrete'. See ?setType.", sep = "" ) 
                    stop( e )
                }                
            )
        } else { 
			.Object <- .ddatamoments(value = value)
        }
		return(.Object)
}

setMethod( "getMean", "datamoments", 
          function( object ) 
          {
              return( object@mean )
          }
)

setMethod( "getVar", "datamoments", 
          function( object ) 
          {
              return( object@var )
          }
)

setMethod( "getFdata", "datamoments", 
          function( object ) 
          {
              return( object@fdata )
          }
)

