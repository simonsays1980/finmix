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

.mcmcpermfix <- setClass("mcmcpermfix", 
                         representation(
                                        Mperm       = "integer",
                                        parperm     = "list",
                                        logperm     = "list"),
                         validity = function(object) 
                         {
                             ## else: OK
                             TRUE
                         },
                         prototype(Mperm    = integer(),
                                   parperm  = list(),
                                   logperm  = list()
                                   )
)

## Getters ##
setMethod("getMperm", "mcmcpermfix", 
          function(object) 
          {
              return(object@Mperm)
          }
)

setMethod("getParperm", "mcmcpermfix", 
          function(object) 
          {
              return(object@parperm)
          }
)

setMethod("getLogperm", "mcmcpermfix", 
          function(object) 
          {
              return(object@logperm)
          }
)

## No setters as users are not intended to modify these ##
## obect.                                               ##
