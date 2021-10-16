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

#' Finmix `csdatamoments` class
#' 
#' Stores moments for indicators of continuous data. Inherited directly from 
#' the [sdatamoments-class] class. 
#' 
#' @slot B A vector storing the between-group heterogeneity.
#' @slot W A vector storing the within-group heterogeneity.
#' @slot T A vector storing the total variance.
#' @slot R A numeric storing the coefficient of determination for univariate 
#'   data.
#' @slot Rdet A numeric storing the coefficient of determination using the 
#'   trace for multivariate data.
#' @slot Rtr A numeric storing the coefficient of determination using the 
#'   determinants for multivariate data. 
#' @exportClass csdatamoments
#' @rdname csdatamoments-class
#' @keywords internal
#' @seealso 
#' * [datamoments-class] for the base class for data moments
#' * [datamoments()] for the constructor of any object of the `datamoments` 
#'   class family
#' * [csdatamoments-class] for the corresponding class defining
#'   moments for data from a discrete-valued finite mixture
.csdatamoments <- setClass("csdatamoments",
  representation(
    B = "vector",
    W = "vector",
    T = "vector",
    R = "numeric",
    Rtr = "numeric",
    Rdet = "numeric"
  ),
  contains = c("sdatamoments"),
  validity = function(object) {
    ## else: ok
    TRUE
  },
  prototype(
    B = vector("numeric"),
    W = vector("numeric"),
    T = vector("numeric"),
    R = numeric(),
    Rtr = numeric(),
    Rdet = numeric()
  )
)

#' Finmix class union of `csdatamoments` and `NULL`
#' 
#' @description
#' Defines a class union such that the object held by a child class can also
#' be `NULL`.
#' 
#' @exportClass csdatamomentsOrNULL
#' @noRd
setClassUnion("csdatamomentsOrNULL", members = c("csdatamoments", "NULL"))

#' Initializer of the `csdatamoments` class
#' 
#' @description
#' Only used implicitly. The initializer calls a function `generateMoments()`
#' to generate in the initialization step the moments for a passed-in `fdata` 
#' object.
#' 
#' @param .Object An object: see the "initialize Methods" section in 
#'   [initialize].
#' @param ... Arguments to specify properties of the new object, to be passed 
#'   to `initialize()`.
#' @param model A finmix [fdata-class] object containing the observations.
#' @noRd
#' 
#' @seealso 
#' * [Classes_Details] for details of class definitions, and 
#' * [setOldClass] for the relation to S3 classes
setMethod(
  "initialize", "csdatamoments",
  function(.Object, ..., value = fdata()) {
    .Object <- callNextMethod(.Object, ..., value = value)
    if (hasY(value) && hasS(value)) {
      .Object <- generateMoments(.Object)
    }
    return(.Object)
  }
)

#' Generate moments for indicators from a mixture with continuous data
#' 
#' @description 
#' Implicit method. Calling [generateMoments()] generates the moments of a
#' finite mixture with continuous data.
#' 
#' @param object An `csdatamoments` object. 
#' @return An `csdatamoments` object with calculated moments.
#' @noRd
setMethod(
  "generateMoments", "csdatamoments",
  function(object) {
    .generateCsdatamoments(object)
  }
)

#' Shows a summary of an `csdatamoments` object.
#' 
#' Calling [show()] on an `csdatamoments` object gives an overview 
#' of the moments of a finite mixture with continuous data.
#' 
#' @param object An `csdatamoments` object.
#' @returns A console output listing the slots and summary information about
#'   each of them. 
#' @exportMethod show
#' @noRd
setMethod(
  "show", "csdatamoments",
  function(object) {
    cat("Object 'sdatamoments'\n")
    cat(
      "     B           : Vector of",
      length(object@B), "\n"
    )
    cat(
      "     W           : Vector of",
      length(object@W), "\n"
    )
    cat(
      "     T           : Vector of",
      length(object@T), "\n"
    )
    if (object@fdata@r > 1) {
      cat("     Rdet        :", object@Rdet, "\n")
      cat("     Rtr         :", object@Rtr, "\n")
    }
    cat(
      "     gmoments    : Object of class",
      class(object@gmoments), "\n"
    )
    cat(
      "     fdata        : Object of class",
      class(object@fdata), "\n"
    )
  }
)

## Getters ##
#' Getter method of `csdatamoments` class.
#' 
#' Returns the `gmoments` slot.
#' 
#' @param object An `csdatamoments` object.
#' @returns The `gmoments` slot of the `object`.
#' @exportMethod getGmoments
#' @noRd
#' 
#' @examples 
#' # Generate an exponential mixture model with two components.
#' f_model <- model("exponential", par = list(lambda = c(0.3, 0.7)), K = 2)
#' # Simulate data from the model.
#' f_data <- simulate(f_model)
#' # Calculate the mixture moments.
#' f_sdatamoms <- sdatamoments(f_data)
#' # Get the moments for the included indicators of the data. 
#' getGmoments(f_sdatamoms)
#' 
#' @seealso 
#' * [datamoments-class] for the base class for model moments
#' * [datamoments()] for the constructor of the `datamoments` 
#'   class family
#' * [csdatamoments-class] for the class definition
#' * [sdatamoments()] for the constructor of the class
setMethod(
  "getGmoments", "csdatamoments",
  function(object) {
    return(object@gmoments)
  }
)

#' Getter method of `csdatamoments` class.
#' 
#' Returns the `WK` slot.
#' 
#' @param object An `csdatamoments` object.
#' @returns The `WK` slot of the `object`.
#' @exportMethod getWK
#' @noRd
#' 
#' @examples 
#' # Generate an exponential mixture model with two components.
#' f_model <- model("exponential", par = list(lambda = c(0.3, 0.7)), K = 2)
#' # Simulate data from the model.
#' f_data <- simulate(f_model)
#' # Calculate the mixture moments.
#' f_sdatamoms <- sdatamoments(f_data)
#' # Get the moments for the included indicators of the data. 
#' getWK(f_sdatamoms)
#' 
#' @seealso 
#' * [datamoments-class] for the base class for model moments
#' * [datamoments()] for the constructor of the `datamoments` 
#'   class family
#' * [csdatamoments-class] for the class definition
#' * [sdatamoments()] for the constructor of the class
setMethod(
  "getWK", "csdatamoments",
  function(object) {
    return(object@WK)
  }
)

#' Getter method of `csdatamoments` class.
#' 
#' Returns the `var` slot.
#' 
#' @param object An `csdatamoments` object.
#' @returns The `var` slot of the `object`.
#' @exportMethod getVar
#' @noRd
#' 
#' @examples 
#' # Generate an exponential mixture model with two components.
#' f_model <- model("exponential", par = list(lambda = c(0.3, 0.7)), K = 2)
#' # Simulate data from the model.
#' f_data <- simulate(f_model)
#' # Calculate the mixture moments.
#' f_sdatamoms <- sdatamoments(f_data)
#' # Get the moments for the included indicators of the data. 
#' getVar(f_sdatamoms)
#' 
#' @seealso 
#' * [datamoments-class] for the base class for model moments
#' * [datamoments()] for the constructor of the `datamoments` 
#'   class family
#' * [csdatamoments-class] for the class definition
#' * [sdatamoments()] for the constructor of the class
setMethod(
  "getVar", "csdatamoments",
  function(object) {
    return(object@var)
  }
)

#' Getter method of `csdatamoments` class.
#' 
#' Returns the `B` slot.
#' 
#' @param object An `csdatamoments` object.
#' @returns The `B` slot of the `object`.
#' @exportMethod getB
#' 
#' @examples 
#' # Generate an exponential mixture model with two components.
#' f_model <- model("exponential", par = list(lambda = c(0.3, 0.7)), K = 2)
#' # Simulate data from the model.
#' f_data <- simulate(f_model)
#' # Calculate the mixture moments.
#' f_sdatamoms <- sdatamoments(f_data)
#' # Get the moments for the included indicators of the data. 
#' getB(f_sdatamoms)
#' 
#' @seealso 
#' * [datamoments-class] for the base class for model moments
#' * [datamoments()] for the constructor of the `datamoments` 
#'   class family
#' * [csdatamoments-class] for the class definition
#' * [sdatamoments()] for the constructor of the class
setMethod(
  "getB", "csdatamoments",
  function(object) {
    return(object@B)
  }
)

#' Getter method of `csdatamoments` class.
#' 
#' Returns the `W` slot.
#' 
#' @param object An `csdatamoments` object.
#' @returns The `W` slot of the `object`.
#' @exportMethod getW
#' @noRd
#' 
#' @examples 
#' # Generate an exponential mixture model with two components.
#' f_model <- model("exponential", par = list(lambda = c(0.3, 0.7)), K = 2)
#' # Simulate data from the model.
#' f_data <- simulate(f_model)
#' # Calculate the mixture moments.
#' f_sdatamoms <- sdatamoments(f_data)
#' # Get the moments for the included indicators of the data. 
#' getW(f_sdatamoms)
#' 
#' @seealso 
#' * [datamoments-class] for the base class for model moments
#' * [datamoments()] for the constructor of the `datamoments` 
#'   class family
#' * [csdatamoments-class] for the class definition
#' * [sdatamoments()] for the constructor of the class
setMethod(
  "getW", "csdatamoments",
  function(object) {
    return(object@W)
  }
)

#' Getter method of `csdatamoments` class.
#' 
#' Returns the `T` slot.
#' 
#' @param object An `csdatamoments` object.
#' @returns The `T` slot of the `object`.
#' @exportMethod getT
#' @noRd
#' 
#' @examples 
#' # Generate an exponential mixture model with two components.
#' f_model <- model("exponential", par = list(lambda = c(0.3, 0.7)), K = 2)
#' # Simulate data from the model.
#' f_data <- simulate(f_model)
#' # Calculate the mixture moments.
#' f_sdatamoms <- sdatamoments(f_data)
#' # Get the moments for the included indicators of the data. 
#' getT(f_sdatamoms)
#' 
#' @seealso 
#' * [datamoments-class] for the base class for model moments
#' * [datamoments()] for the constructor of the `datamoments` 
#'   class family
#' * [csdatamoments-class] for the class definition
#' * [sdatamoments()] for the constructor of the class
setMethod(
  "getT", "csdatamoments",
  function(object) {
    return(object@T)
  }
)

#' Getter method of `csdatamoments` class.
#' 
#' Returns the `R` slot.
#' 
#' @param object An `csdatamoments` object.
#' @returns The `R` slot of the `object`.
#' @exportMethod getR
#' @noRd
#' 
#' @examples 
#' # Generate an exponential mixture model with two components.
#' f_model <- model("exponential", par = list(lambda = c(0.3, 0.7)), K = 2)
#' # Simulate data from the model.
#' f_data <- simulate(f_model)
#' # Calculate the mixture moments.
#' f_sdatamoms <- sdatamoments(f_data)
#' # Get the moments for the included indicators of the data. 
#' getR(f_sdatamoms)
#' 
#' @seealso 
#' * [datamoments-class] for the base class for model moments
#' * [datamoments()] for the constructor of the `datamoments` 
#'   class family
#' * [csdatamoments-class] for the class definition
#' * [sdatamoments()] for the constructor of the class
setMethod(
  "getR", "csdatamoments",
  function(object) {
    return(object@R)
  }
)

#' Getter method of `csdatamoments` class.
#' 
#' Returns the `Rtr` slot.
#' 
#' @param object An `csdatamoments` object.
#' @returns The `Rtr` slot of the `object`.
#' @exportMethod getRtr
#' @noRd
#' 
#' @examples 
#' # Generate an exponential mixture model with two components.
#' f_model <- model("exponential", par = list(lambda = c(0.3, 0.7)), K = 2)
#' # Simulate data from the model.
#' f_data <- simulate(f_model)
#' # Calculate the mixture moments.
#' f_sdatamoms <- sdatamoments(f_data)
#' # Get the moments for the included indicators of the data. 
#' getRtr(f_sdatamoms)
#' 
#' @seealso 
#' * [datamoments-class] for the base class for model moments
#' * [datamoments()] for the constructor of the `datamoments` 
#'   class family
#' * [csdatamoments-class] for the class definition
#' * [sdatamoments()] for the constructor of the class
setMethod(
  "getRtr", "csdatamoments",
  function(object) {
    return(object@Rtr)
  }
)

#' Getter method of `csdatamoments` class.
#' 
#' Returns the `Rdet` slot.
#' 
#' @param object An `csdatamoments` object.
#' @returns The `Rdet` slot of the `object`.
#' @exportMethod getRdet
#' @noRd
#' 
#' @examples 
#' # Generate an exponential mixture model with two components.
#' f_model <- model("exponential", par = list(lambda = c(0.3, 0.7)), K = 2)
#' # Simulate data from the model.
#' f_data <- simulate(f_model)
#' # Calculate the mixture moments.
#' f_sdatamoms <- sdatamoments(f_data)
#' # Get the moments for the included indicators of the data. 
#' getRdet(f_sdatamoms)
#' 
#' @seealso 
#' * [datamoments-class] for the base class for model moments
#' * [datamoments()] for the constructor of the `datamoments` 
#'   class family
#' * [csdatamoments-class] for the class definition
#' * [sdatamoments()] for the constructor of the class
setMethod(
  "getRdet", "csdatamoments",
  function(object) {
    return(object@Rdet)
  }
)

#' Getter method of `csdatamoments` class.
#' 
#' Returns the `fdata` slot.
#' 
#' @param object An `csdatamoments` object.
#' @returns The `fdata` slot of the `object`.
#' @exportMethod getFdata
#' @noRd
#' 
#' @examples 
#' # Generate an exponential mixture model with two components.
#' f_model <- model("exponential", par = list(lambda = c(0.3, 0.7)), K = 2)
#' # Simulate data from the model.
#' f_data <- simulate(f_model)
#' # Calculate the mixture moments.
#' f_sdatamoms <- sdatamoments(f_data)
#' # Get the moments for the included indicators of the data. 
#' getFdata(f_sdatamoms)
#' 
#' @seealso 
#' * [datamoments-class] for the base class for model moments
#' * [datamoments()] for the constructor of the `datamoments` 
#'   class family
#' * [csdatamoments-class] for the class definition
#' * [sdatamoments()] for the constructor of the class
setMethod(
  "getFdata", "csdatamoments",
  function(object) {
    return(object@fdata)
  }
)

## Setters ##
## No setters, as it users are adviced not to manipulate moment objects ##

### Private functions
### These functions are not exported
#' Generate data moments for finite mixture data
#' 
#' @description 
#' Only called implicitly. generates all moments of finite mixture data in a 
#' `fdata` object.
#' 
#' @param object A `csdatamoments` object to contain all calculated
#'   moments. 
#' @returns A `csdatamoments` object containing all moments of the 
#'   inite mixture data.
#' @noRd
".generateCsdatamoments" <- function(object) {
  ## enforce column.wise ordering ##
  datam <- getColY(object@fdata)
  classm <- getColS(object@fdata)
  ## Calculate the between-group variance ##
  ## 'B' is an r x r matrix ##
  gmeans <- object@gmoments@mean
  nkm <- object@gmoments@NK
  ## Calculate the total heterogeneity ##
  ## 'T' is an r x r array ##
  object@T <- var(datam, na.rm = TRUE) * nrow(datam)
  ## Calculate the within-group heterogeneity ##
  ## 'W' is an r x r array ##
  wkm <- object@gmoments@WK
  object@W <- apply(wkm, c(1, 2), sum, na.rm = TRUE)
  ## Calculate between-group heterogeneity ##
  ## 'B' is an r x r array ##
  object@B <- object@T - object@W
  ## Calculate coefficient of determination ##
  ## 'Rtr' is an 1 x 1 numeric ##
  ## 'Rdet' is an 1 x 1 numeric ##
  if (object@fdata@r > 1) {
    r <- NA
    object@R <- as.numeric(r)
    object@Rtr <- 1 - sum(diag(object@W), na.rm = TRUE) /
      sum(diag(object@T), na.rm = TRUE)
    object@Rdet <- 1 - det(object@W) / det(object@T)
  } else {
    rtr <- NA
    rdet <- NA
    object@Rtr <- as.numeric(rtr)
    object@Rdet <- as.numeric(rdet)
    object@R <- 1 - object@W[1] / object@T[1]
  }
  return(object)
}
