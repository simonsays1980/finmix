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
#' Finmix `mcmc` class
#' 
#' @description 
#' This class defines hyper-parameters for the MCMC procedure. This is a main 
#' class of the `finmix` package that must be defined for estimating a finite 
#' mixture model. 
#' 
#' @slot burnin An integer defining the number of steps in the burn-in phase of 
#'   Gibbs-sampling.
#' @slot M An integer defining the number of steps in Gibbs-sampling to be 
#'   stored. 
#' @slot startpar A logical indicating, if starting by sampling the 
#'   parameters. If `FALSE` sampling starts by sampling the indicators `S`.
#' @slot storeS An integer specifying how many of the last sampled indicators 
#'   should be stored in the output. 
#' @slot storepost A logical indicating if the posterior probabilities should 
#'   be stored. This becomes for example important for specific relabeling 
#'   algorithms, but also for analysis.
#' @slot ranperm A logical indicating, if random permutation should be used. If 
#'   `TRUE` the parameters are permutated randomly between the number of 
#'   components after each sampling step in MCMC.
#' @slot storeinv A logical indicating if the inverse variance-covariance 
#'   matrices for multivariate normal or Student-t mixtures should be stored. 
#' @exportClass mcmc
#' @rdname mcmc-class
#' 
#' @seealso 
#' * [mcmc()] for the class constructor
#' * [mcmcstart()] for completion of slots 
#' * [mixturemcmc()] for further information about the MCMC sampling
.mcmc <- setClass("mcmc",
  representation(
    burnin = "integer",
    M = "integer",
    startpar = "logical",
    storeS = "integer",
    storepost = "logical",
    ranperm = "logical",
    storeinv = "logical"
  ),
  validity = function(object) {
    .valid.MCMC(object)
    ## else: OK
    TRUE
  },
  prototype(
    burnin = integer(),
    M = integer(),
    startpar = logical(),
    storeS = integer(),
    storepost = logical(),
    ranperm = logical(),
    storeinv = logical()
  )
)

#' Constructor for `mcmc` class
#' 
#' @description 
#' Calling [mcmc()] constructs an object of class `mcmc` that specifies the 
#' hyper-parameters for the MCMC procedure. Each MCMC sampling needs an `mcmc` 
#' object that specifies the way, how MCMC sampling should be performed and 
#' what kind and how much of data should be stored.
#' 
#' @param burnin An integer defining the number of steps in the burn-in phase of 
#'   Gibbs-sampling.
#' @param M An integer defining the number of steps in Gibbs-sampling to be 
#'   stored. 
#' @param startpar A logical indicating, if starting by sampling the 
#'   parameters. If `FALSE` sampling starts by sampling the indicators `S`.
#' @param storeS An integer specifying how many of the last sampled indicators 
#'   should be stored in the output. 
#' @param storepost A logical indicating if the posterior probabilities should 
#'   be stored. This becomes for example important for specific relabeling 
#'   algorithms, but also for analysis.
#' @param ranperm A logical indicating, if random permutation should be used. If 
#'   `TRUE` the parameters are permutated randomly between the number of 
#'   components after each sampling step in MCMC.
#' @param storeinv A logical indicating if the inverse variance-covariance 
#'   matrices for multivariate normal or Student-t mixtures should be stored.
#' @return An object of class `mcmc` containing all hyper-parameters for MCMC
#'   sampling.
#' @export
#' @name mcmc
#' 
#' @examples 
#' f_mcmc <- mcmc()
#' 
#' @seealso 
#' * [mcmc-class] for the definition of the `mcmc` class
#' * [mcmcstart()] for setting up all objects for MCMC sampling
#' * [mixturemcmc()] for running MCMC sampling for finite mixture models
"mcmc" <- function(burnin = 0, M = 5000,
                   startpar = TRUE, storeS = 1000,
                   storepost = TRUE, ranperm = TRUE,
                   storeinv = TRUE) {
  .mcmc(
    burnin = as.integer(burnin),
    M = as.integer(M), startpar = startpar,
    storeS = as.integer(storeS), storepost = storepost,
    ranperm = ranperm, storeinv = storeinv
  )
}

#' Shows a summary of an `mcmc` object.
#' 
#' Calling [show()] on an `mcmc` object gives an overview 
#' of the `mcmc` object.
#' 
#' @param object A `mcmc` object.
#' @returns A console output listing the slots and summary information about
#'   each of them. 
#' @exportMethod show
#' @noRd
#' 
#' @seealso 
#' * [mcmc-class] for the class definition
#' * [mcmc()] for the constructor of the class 
setMethod(
  "show", "mcmc",
  function(object) {
    cat("Object 'mcmc'\n")
    cat("     class       :", class(object), "\n")
    cat("     burnin      :", object@burnin, "\n")
    cat("     M           :", object@M, "\n")
    cat("     startpar    :", object@startpar, "\n")
    cat("     storeS      :", object@storeS, "\n")
    cat("     storepost   :", object@storepost, "\n")
    cat("     ranperm     :", object@ranperm, "\n")
    cat("     storeinv    :", object@storeinv, "\n")
  }
)

## Getters ##
#' Getter method of `mcmc` class.
#' 
#' Returns the `burnin` slot.
#' 
#' @param object An `mcmc` object.
#' @returns The `burnin` slot of the `object`.
#' @noRd
#' @export
#' @examples 
#' # Generate an mcmc object
#' f_mcmc <- mcmc()
#' # Get the slot
#' getBurnin(f_mcmc)
#' 
#' @seealso 
#' * [mcmc-class] for the class definition
#' * [mcmc()] for the constructor of the `mcmc` class
setMethod(
  "getBurnin", "mcmc",
  function(object) {
    return(object@burnin)
  }
)

#' Getter method of `mcmc` class.
#' 
#' Returns the `M` slot.
#' 
#' @param object An `mcmc` object.
#' @returns The `M` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' # Generate an mcmc object
#' f_mcmc <- mcmc()
#' # Get the slot
#' getM(f_mcmc)
#' 
#' @seealso 
#' * [mcmc-class] for the class definition
#' * [mcmc()] for the constructor of the `mcmc` class
setMethod(
  "getM", "mcmc",
  function(object) {
    return(object@M)
  }
)

#' Getter method of `mcmc` class.
#' 
#' Returns the `startpar` slot.
#' 
#' @param object An `mcmc` object.
#' @returns The `startpar` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' # Generate an mcmc object
#' f_mcmc <- mcmc()
#' # Get the slot
#' getStartpar(f_mcmc)
#' 
#' @seealso 
#' * [mcmc-class] for the class definition
#' * [mcmc()] for the constructor of the `mcmc` class
setMethod(
  "getStartpar", "mcmc",
  function(object) {
    return(object@startpar)
  }
)

#' Getter method of `mcmc` class.
#' 
#' Returns the `storeS` slot.
#' 
#' @param object An `mcmc` object.
#' @returns The `storeS` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' # Generate an mcmc object
#' f_mcmc <- mcmc()
#' # Get the slot
#' getStoreS(f_mcmc)
#' 
#' @seealso 
#' * [mcmc-class] for the class definition
#' * [mcmc()] for the constructor of the `mcmc` class
setMethod(
  "getStoreS", "mcmc",
  function(object) {
    return(object@storeS)
  }
)

#' Getter method of `mcmc` class.
#' 
#' Returns the `storepost` slot.
#' 
#' @param object An `mcmc` object.
#' @returns The `storepost` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' # Generate an mcmc object
#' f_mcmc <- mcmc()
#' # Get the slot
#' getStorepost(f_mcmc)
#' 
#' @seealso 
#' * [mcmc-class] for the class definition
#' * [mcmc()] for the constructor of the `mcmc` class
setMethod(
  "getStorepost", "mcmc",
  function(object) {
    return(object@storepost)
  }
)

#' Getter method of `mcmc` class.
#' 
#' Returns the `ranperm` slot.
#' 
#' @param object An `mcmc` object.
#' @returns The `ranperm` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' # Generate an mcmc object
#' f_mcmc <- mcmc()
#' # Get the slot
#' getRanperm(f_mcmc)
#' 
#' @seealso 
#' * [mcmc-class] for the class definition
#' * [mcmc()] for the constructor of the `mcmc` class
setMethod(
  "getRanperm", "mcmc",
  function(object) {
    return(object@ranperm)
  }
)

## Setters ##
#' Setter method of `mcmc` class.
#' 
#' Sets the slot. Returns the none.
#' 
#' @param object An `mcmc` object.
#' @param value An integer defining the new value for the `@@burnin` slot.
#' @returns None.
#' @noRd
#' 
#' @examples 
#' # Generate an mcmc object
#' f_mcmc <- mcmc()
#' # Set the slot
#' setBurnin(f_mcmc) <- as.integer(2000)
#' 
#' @seealso 
#' * [mcmc-class] for the class definition
#' * [mcmc()] for the constructor of the `mcmc` class
setReplaceMethod(
  "setBurnin", "mcmc",
  function(object, value) {
    object@burnin <- as.integer(value)
    validObject(object)
    return(object)
  }
)

#' Setter method of `mcmc` class.
#' 
#' Sets the slot. Returns the none.
#' 
#' @param object An `mcmc` object.
#' @param value An integer defining the new value for the `@@M` slot.
#' @returns None.
#' @noRd
#' 
#' @examples 
#' # Generate an mcmc object
#' f_mcmc <- mcmc()
#' # Set the slot
#' setM(f_mcmc) <- as.integer(20000)
#' 
#' @seealso 
#' * [mcmc-class] for the class definition
#' * [mcmc()] for the constructor of the `mcmc` class
setReplaceMethod(
  "setM", "mcmc",
  function(object, value) {
    object@M <- as.integer(value)
    validObject(object)
    return(object)
  }
)

#' Setter method of `mcmc` class.
#' 
#' Sets the slot. Returns the none.
#' 
#' @param object An `mcmc` object.
#' @param value An integer defining the new value for the `@@startpar` slot.
#' @returns None.
#' @noRd
#' 
#' @examples 
#' # Generate an mcmc object
#' f_mcmc <- mcmc()
#' # Set the slot
#' setStartpar(f_mcmc) <- FALSE
#' 
#' @seealso 
#' * [mcmc-class] for the class definition
#' * [mcmc()] for the constructor of the `mcmc` class
setReplaceMethod(
  "setStartpar", "mcmc",
  function(object, value) {
    object@startpar <- value
    validObject(object)
    return(object)
  }
)

#' Setter method of `mcmc` class.
#' 
#' Sets the slot. Returns the none.
#' 
#' @param object An `mcmc` object.
#' @param value An integer defining the new value for the `@@storeS` slot.
#' @returns None.
#' @noRd
#' 
#' @examples 
#' # Generate an mcmc object
#' f_mcmc <- mcmc()
#' # Set the slot
#' setStoreS(f_mcmc) <- as.integer(500)
#' 
#' @seealso 
#' * [mcmc-class] for the class definition
#' * [mcmc()] for the constructor of the `mcmc` class
setReplaceMethod(
  "setStoreS", "mcmc",
  function(object, value) {
    object@storeS <- as.integer(value)
    validObject(object)
    return(object)
  }
)

#' Setter method of `mcmc` class.
#' 
#' Sets the slot. Returns the none.
#' 
#' @param object An `mcmc` object.
#' @param value An integer defining the new value for the `@@storepost` slot.
#' @returns None.
#' @exportMethod setStorepost<-
#' @noRd
#' 
#' @examples 
#' # Generate an mcmc object
#' f_mcmc <- mcmc()
#' # Set the slot
#' setStorepost(f_mcmc) <- FALSE
#' 
#' @seealso 
#' * [mcmc-class] for the class definition
#' * [mcmc()] for the constructor of the `mcmc` class
setReplaceMethod(
  "setStorepost", "mcmc",
  function(object, value) {
    object@storepost <- value
    validObject(object)
    return(object)
  }
)

#' Setter method of `mcmc` class.
#' 
#' Sets the slot. Returns the none.
#' 
#' @param object An `mcmc` object.
#' @param value An integer defining the new value for the `@@ranperm` slot.
#' @returns None.
#' @noRd
#' 
#' @examples 
#' # Generate an mcmc object
#' f_mcmc <- mcmc()
#' # Set the slot
#' setRanperm(f_mcmc) <- FALSE
#' 
#' @seealso 
#' * [mcmc-class] for the class definition
#' * [mcmc()] for the constructor of the `mcmc` class
setReplaceMethod(
  "setRanperm", "mcmc",
  function(object, value) {
    object@ranperm <- value
    validObject(object)
    return(object)
  }
)

### Private functions
### These functions are not exported

### Valid mcmc: The number of burnins @burnin and the number of
### last indicator vectors to store @storeS must be non-negative
### 'integers'. The number of MCMC draws @M must be a positive 'integer'.
#' Check validity of an `mcmc` object
#' 
#' @description 
#' For internal usage only. This function checks if the different slots of the
#' `mcmc` object are valid. It checks if slots `@@burnin`, `@@M`, and 
#' `@@storeS` are set to non-negative values and if slot `@@storeS` does not 
#' call for more indicators to store than iterations in the MCMC sampling.
#' 
#' @param object An `mcmc` object to be checked.
#' @return None. If checks do not pass an error is thrown.
#' @noRd
#' @seealso 
#' * [mcmc()] for the calling function
".valid.MCMC" <- function(object) {
  if (object@burnin < as.integer(0)) {
    stop(paste("Number of Burn-In draws in slot 'burnin' must be ",
      "nonnegative.",
      sep = ""
    ))
  } else if (object@M <= as.integer(0)) {
    stop("Number of MCMC draws in slot 'M' must be positive. ")
  } else if (object@storeS < as.integer(0)) {
    stop(paste("Number of indicators to store in slot 'storeS' must be ",
      "nonnegative.",
      sep = ""
    ))
  }
  if (object@storeS > object@M) {
    stop(paste("Number of indicators to store in slot 'storeS' must be ",
      "smaller or equal to the number of MCMC draws in slot 'M'.",
      sep = ""
    ))
  }
}
