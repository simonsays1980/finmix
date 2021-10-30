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

#' Finmix `prior` class
#' 
#' The `prior` class stores the specifications for the prior distribution used 
#' for Bayesian estimation of the finite mixture parameters and weights. There 
#' exists next to the general constructor also an advanced constructor that 
#' specifies a data dependent prior. See [priordefine()] for this advanced 
#' constructor. 
#' 
#' @slot weight A matrix storing the prior parameters for the `weight` of a 
#'   finite mixture model.
#' @slot par A list storing the prior parameters for the parameters of a finite 
#'   mixture model.
#' @slot type A character specifying what type of prior should be used in 
#'   Bayesian estimation. Either `"independent"` for an independent prior 
#'   distribution or `"condconjugate"` for a conditionally conjugate prior 
#'   distribution.
#' @slot hier A logical defining, if the used prior should be hierarchical. 
#'   Hierarchical prior are often more robust, but need an additional layer in 
#'   sampling, so computing costs increase.
#' @exportClass prior
#' @rdname prior-class
#' 
#' @seealso 
#' * [prior()] for the general constructor of this class
#' * [priordefine()] for the advanced constructor of this class
#' 
#' @references 
#' * Frühwirth-Schnatter, S (2006), "Finite Mixture and Markov Switching Models"
.prior <- setClass("prior",
  representation(
    weight = "matrix",
    par = "list",
    type = "character",
    hier = "logical"
  ),
  validity = function(object) {
    .valid.type.Prior(object)
    ## else: OK
    TRUE
  },
  prototype(
    weight = matrix(),
    par = list(),
    type = character(),
    hier = logical()
  )
)

### ----------------------------------------------------------------
### Constructors
### ----------------------------------------------------------------

#' Constructor for `prior` class
#' 
#' @description 
#' Calling [prior()] constructs an object of class [prior][prior-class]. The 
#' constructor can be called without providing any arguments, but the prior 
#' has to be filled with appropriate parameters when MCMC sampling should be 
#' performed. 
#' 
#' There exists next to the general constructor also an advanced constructor 
#' that specifies a data dependent prior. See [priordefine()] for this advanced 
#' constructor.
#' 
#' @slot weight A matrix storing the prior parameters for the `weight` of a 
#'   finite mixture model.
#' @slot par A list storing the prior parameters for the parameters of a finite 
#'   mixture model.
#' @slot type A character specifying what type of prior should be used in 
#'   Bayesian estimation. Either `"independent"` for an independent prior 
#'   distribution or `"condconjugate"` for a conditionally conjugate prior 
#'   distribution.
#' @slot hier A logical defining, if the used prior should be hierarchical. 
#'   Hierarchical prior are often more robust, but need an additional layer in 
#'   sampling, so computing costs increase.
#' @export
#' @name prior
#' 
#' @examples 
#' # Call the default constructor without any arguments.
#' f_prior <- prior()
#' 
#' @seealso 
#' * [prior()] for the general constructor of this class
#' * [priordefine()] for the advanced constructor of this class
#' 
#' @references 
#' * Fr\"uhwirth-Schnatter, S (2006), "Finite Mixture and Markov Switching Models"
"prior" <- function(weight = matrix(), par = list(),
                    type = c("independent", "condconjugate"),
                    hier = TRUE) {
  type <- match.arg(type)
  .prior(
    weight = weight, par = par,
    type = type, hier = hier
  )
}

#' Advanced constructor for the `prior` class
#' 
#' This constructor defines a data dependent prior with parameters by matching 
#' moments. As a consequence it needs as inputs an `fdata` object and a `model` 
#' object. The prior distributions chosen and the methods how parameters are 
#' computed are described in Frühwirth-Schnatter (2006). 
#' 
#' @param fdata An `fdata` object holding the data. Observations in slot `@@y` 
#'   must be existent.
#' @param model A `model` object specifying the finite mixture model. 
#' @param varargin `NULL` or a `prior` object. This enables the user to pass in 
#'   an already constructed prior object that gets then completed.
#' @param prior.wagner A logical indicating, if the prior from Wagner (2007) 
#'   should be used in case of an exponential mixture model. 
#' @param s A numeric specifying the standard deviation `s` for the 
#'   Metropolis-Hastings proposal.  
#' @return A fully specified `prior` object.
#' @export
#' @name priordefine
#' 
#' @examples 
#' # Create a Poisson mixture model.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 0.7)), K = 2)
#' # Simulate data from the model.
#' f_data <- simulate(f_model)
#' # Use the advanced constructor to generate a prior.
#' f_prior <- priordefine(f_data, f_model)
#' 
#' @seealso 
#' * [prior-class] for the class definition
#' * [prior()] for the default constructor of the class
#' 
#' @references 
#' * Fr\"uwirth-Schnatter, S. (2006), "Finite Mixture and Markov Switching 
#'   Models"
#' * Wagner, H. (2007), "Bayesian analysis of mixtures of exponentials", 
#'   Journal of Applied Mathematics, Statistics and Informatics 3, 165-183
"priordefine" <- function(fdata = fdata(), model = model(),
                          varargin = NULL, prior.wagner = TRUE, s = 5.0) {
  .check.fdata.model.Prior(fdata, model)
  if (!is.null(varargin)) {
    .check.varargin.Prior(varargin)
  }
  object <- .prior(hier = TRUE, type = "independent")
  generatePrior(object,
    fdata = fdata, model = model,
    varargin = varargin, prior.wagner = prior.wagner, s
  )
}

### ==================================================================
### Has methods
### ------------------------------------------------------------------
#' Checks for parameters in a `prior` object
#' 
#' @description 
#' Calling `hasPriorPar()` checks if `model`-appropriate parameters are stored 
#' in the `prior` object.
#' 
#' @param object A `prior` object containing the specifications for the prior.
#' @param model A `model` object containing the specifications for the model.
#' @param verbose A logical indicating, if the output should be verbose.
#' @exportMethod hasPriorPar
#' @keywords internal 
#' 
#' @examples 
#' # Define a Poisson mixture model.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 0.7)), K = 2)
#' # Call the default constructor.
#' f_prior <- prior()
#' # Check if the prior has appropriate parameters defined.
#' hasPriorPar(f_prior, f_model)
#' \dontrun{hasPriorPar(f_prior, f_model, TRUE)}
#' 
#' @seealso 
#' * [prior-class] for the definition of the `prior` class
#' * [model-class] for the definition of the `model` class
setMethod(
  "hasPriorPar", signature(
    object = "prior",
    model = "model",
    verbose = "ANY"
  ),
  function(object, model, verbose = FALSE) {
    .haspar.Prior(object, model, verbose)
  }
)

#' Checks for parameters in a `prior` object
#' 
#' @description 
#' Calling `hasPriorWeight()` checks if `model`-appropriate weight parameters 
#' are stored in the `prior` object.
#' 
#' @param object A `prior` object containing the specifications for the prior.
#' @param model A `model` object containing the specifications for the model.
#' @param verbose A logical indicating, if the output should be verbose.
#' @exportMethod hasPriorWeight
#' @rdname hasPriorWeight
#' @keywords internal
#' 
#' @examples 
#' # Define a Poisson mixture model.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 0.7)), K = 2)
#' # Call the default constructor.
#' f_prior <- prior()
#' # Check if the prior has appropriate parameters defined.
#' hasPriorWeight(f_prior, f_model)
#' \dontrun{hasPriorWeight(f_prior, f_model, TRUE)}
#' 
#' @seealso 
#' * [prior-class] for the definition of the `prior` class
#' * [model-class] for the definition of the `model` class
setMethod(
  "hasPriorWeight", signature(
    object = "prior",
    model = "model",
    verbose = "ANY"
  ),
  function(object, model, verbose = FALSE) {
    if (!all(is.na(object@weight))) {
      if (ncol(object@weight) == model@K) {
        return(TRUE)
      } else {
        if (verbose) {
          stop(paste("Wrong dimension of ",
            "slot 'weight' of ",
            "'prior' object. ",
            "Weights must be of ",
            "dimension 1 x K.",
            sep = ""
          ))
        } else {
          return(FALSE)
        }
      }
    } else {
      if (verbose) {
        stop(paste("Slot 'weight' of 'prior' ",
          "object is empty.",
          sep = ""
        ))
      } else {
        return(FALSE)
      }
    }
  }
)

#' Generates `prior` object
#' 
#' @description 
#' Calling `generatePrior()` generates the `prior` object when [priordefine()] 
#' had been called. When this function is called all checks have been passed 
#' and `prior` construction can take place. 
#' 
#' @param object A `prior` object to store the prior parameters and weights. 
#' @param fdata An `fdata` object holding the data. Observations in slot `@@y` 
#'   must be existent.
#' @param model A `model` object specifying the finite mixture model. 
#' @param varargin `NULL` or a `prior` object. This enables the user to pass in 
#'   an already constructed prior object that gets then completed.
#' @param prior.wagner A logical indicating, if the prior from Wagner (2007) 
#'   should be used in case of an exponential mixture model. 
#' @param s A numeric specifying the standard deviation `s` for the 
#'   Metropolis-Hastings proposal. 
#' @rdname generatePrior
#' @noRd
#' 
#' @seealso 
#' * [prior-class] for the class definition
#' * [priordefine()] for the advanced class constructor using this method 
setMethod(
  "generatePrior", "prior",
  function(object, fdata, model, varargin, prior.wagner, s) {
    dist <- model@dist
    if (dist == "poisson") {
      object <- .generatePriorPoisson(object, fdata, model,
        varargin = varargin
      )
    } else if (dist == "cond.poisson") {
      object <- .generatePriorCondPoisson(
        object, fdata, model,
        s
      )
    } else if (dist == "binomial") {
      object <- .generatePriorBinomial(object, model)
    } else if (dist == "exponential") {
      object <- .generatePriorExponential(object, fdata, model,
        varargin = varargin,
        prior.wagner
      )
    } else {
      object <- .generatePriorNorstud(object, fdata, model,
        varargin = varargin
      )
      if (dist == "student" || dist == "studmult") {
        object <- .generateDfPrior(object)
        warning(paste("A 'prior' object for a Student-t model ",
          "needs a tuning vector named 'mhtune' ",
          "to be added by the user to slot @par$df ",
          "of the 'prior' object.",
          sep = ""
        ),
        call. = FALSE
        )
      }
    }
    .generatePriorWeight(object, model)
  }
)

#' Shows a summary of a `prior` object.
#' 
#' Calling [show()] on a `prior` object gives an overview 
#' of the slots of a `prior` object. 
#' 
#' @param object A `prior` object.
#' @returns A console output listing the slots and summary information about
#'   each of them. 
#' @exportMethod show
#' @keywords internal
#' @seealso 
#' * [prior-class] for the class definition
#' * [prior()] for the basic constructor of the class
#' * [priordefine()] for the advanced constructor of the class
setMethod(
  "show", "prior",
  function(object) {
    cat("Object 'prior'\n")
    cat("     class       :", class(object), "\n")
    cat("     hier        :", object@hier, "\n")
    cat("     type        :", object@type, "\n")
    cat(
      "     par         : List of",
      length(object@par), "\n"
    )
    if (!all(is.na(object@weight))) {
      cat(
        "     weight      :",
        paste(dim(object@weight), collapse = "x"), "\n"
      )
    }
  }
)

## Getters ##
#' Getter method of `prior` class.
#' 
#' Returns the `weight` slot.
#' 
#' @param object An `prior` object.
#' @returns The `weight` slot of the `object`.
#' @exportMethod getWeight
#' @keywords internal
#' 
#' @examples 
#' # Generate a prior object. 
#' f_prior <- prior()
#' # Get the slot.
#' getWeight(f_prior)
setMethod(
  "getWeight", "prior",
  function(object) {
    return(object@weight)
  }
)

## Getters ##
#' Getter method of `prior` class.
#' 
#' Returns the `par` slot.
#' 
#' @param object An `prior` object.
#' @returns The `par` slot of the `object`.
#' @exportMethod getPar
#' @keywords internal
#' 
#' @examples 
#' # Generate a prior object. 
#' f_prior <- prior()
#' # Get the slot.
#' getPar(f_prior)
setMethod(
  "getPar", "prior",
  function(object) {
    return(object@par)
  }
)

## Getters ##
#' Getter method of `prior` class.
#' 
#' Returns the `type` slot.
#' 
#' @param object An `prior` object.
#' @returns The `type` slot of the `object`.
#' @exportMethod getType
#' @keywords internal
#' 
#' @examples 
#' # Generate a prior object. 
#' f_prior <- prior()
#' # Get the slot.
#' getType(f_prior)
setMethod(
  "getType", "prior",
  function(object) {
    return(object@type)
  }
)

## Getters ##
#' Getter method of `prior` class.
#' 
#' Returns the `hier` slot.
#' 
#' @param object An `prior` object.
#' @returns The `hier` slot of the `object`.
#' @exportMethod getHier
#' @keywords internal
#' 
#' @examples 
#' # Generate a prior object. 
#' f_prior <- prior()
#' # Get the slot.
#' getHier(f_prior)
setMethod(
  "getHier", "prior",
  function(object) {
    return(object@hier)
  }
)

## Setters ##
#' Setter method of `prior` class.
#' 
#' Sets the slot. Returns the none.
#' 
#' @param object An `prior` object.
#' @param value An integer defining the new value for the `@@weight` slot.
#' @returns None.
#' @exportMethod setWeight<-
#' @keywords internal
#' 
#' @examples 
#' # Generate a prior object.
#' f_prior <- prior()
#' # Set the slot.
#' setWeight(f_prior) <- matrix(c(0.5, 0.5), nrow = 1)
setReplaceMethod(
  "setWeight", "prior",
  function(object, value) {
    object@weight <- value
    validObject(object)
    return(object)
  }
)

#' Setter method of `prior` class.
#' 
#' Sets the slot. Returns the none.
#' 
#' @param object An `prior` object.
#' @param value An integer defining the new value for the `@@par` slot.
#' @returns None.
#' @exportMethod setPar<-
#' @keywords internal
#' 
#' @examples 
#' # Generate a prior object.
#' f_prior <- prior()
#' # Set the slot.
#' setPar(f_prior) <- setPar(f_prior) <- list(a = matrix(c(1.2, 0.8), nrow = 1), 
#'                                            b = matrix(c(2.3, 0.4), nrow = 1))
setReplaceMethod(
  "setPar", "prior",
  function(object, value) {
    object@par <- value
    validObject(object)
    return(object)
  }
)

#' Setter method of `prior` class.
#' 
#' Sets the slot. Returns none.
#' 
#' @param object An `prior` object.
#' @param value An integer defining the new value for the `@@type` slot.
#' @returns None.
#' @exportMethod setType<-
#' @keywords internal
#' 
#' @examples 
#' # Generate a prior object.
#' f_prior <- prior()
#' # Set the slot.
#' setType(f_prior) <- "condconjugate"
setReplaceMethod(
  "setType", "prior",
  function(object, value) {
    object@type <- value
    validObject(object)
    return(object)
  }
)

#' Setter method of `prior` class.
#' 
#' Sets the slot. Returns none.
#' 
#' @param object An `prior` object.
#' @param value An integer defining the new value for the `@@hier` slot.
#' @returns None.
#' @exportMethod setHier<-
#' @keywords internal
#' 
#' @examples 
#' # Generate a prior object.
#' f_prior <- prior()
#' # Set the slot.
#' setHier(f_prior) <- TRUE
setReplaceMethod(
  "setHier", "prior",
  function(object, value) {
    object@hier <- value
    validObject(object)
    return(object)
  }
)

### Private functions
### These functions are not exported

#' Check validity of `fdata` and `model` objects for prior generation
#' 
#' @description 
#' For internal usage only. This function checks the validity of the passed in 
#' objects `fdata` and `model` in [priordefine()]. This includes checking, if 
#' the `fdata` object contains observations in slot `@@y` and, if the `model` 
#' object is valid. Finally, the consistency between the two objects is checked. 
#' 
#' @param fdata.obj An `fdata` object. Must contain observations in slot `@@y` 
#'   to pass the checks. 
#' @param model.obj A `model` object. Must be specified by parameters in slot 
#'   `@@par` and number of components in `@@K`.
#' @returns None. If the checks do not pass, an error is thrown.
#' @noRd
".check.fdata.model.Prior" <- function(fdata.obj, model.obj) {
  .valid.Fdata(fdata.obj)
  hasY(fdata.obj, verbose = TRUE)
  .init.valid.Model(model.obj)
  .valid.fdata.model.Prior(fdata.obj, model.obj)
}

#' Check validity of `varargin` object for prior generation
#' 
#' @description 
#' For internal usage only. This function checks the optional `prior` object 
#' passed in to [priordefine()]. This object has to be of class 
#' [model-class] and has to be valid as this. 
#' 
#' @param obj Any object. 
#' @returns None. If the checks do not pass, an error is thrown.
#' @noRd
".check.varargin.Prior" <- function(obj) {
  if (!inherits(obj, "prior")) {
    stop(paste("Argument 'varargin' is not of class 'prior'. ",
      "If argument 'varargin' in 'priordefine()' is ",
      "specified, it must be of class 'prior'.",
      sep = ""
    ))
  } else {
    validObject(obj)
  }
}

### Has
### hasPar Prior
#' Checks for parameters in `prior` object
#' 
#' @description 
#' For internal usage only. This function checks, if a given 
#' [model-class] contains specified parameters in its slot `@@par`.
#' 
#' @param obj A `prior` object to be checked.
#' @param model.obj A `model` object providing the model distribution for 
#'   which prior parameters should be checked. 
#' @param verbose A logical indicating, if the output should be verbose or 
#'   silent.
#' @returns Either a logical indicating, if the passed-in `prior` object 
#'   contains defined parameters or verbose output. Throws an error, if the 
#'   checks do not pass. 
#' @noRd
".haspar.Prior" <- function(obj, model.obj, verbose) {
  dist <- model.obj@dist
  if (dist == "poisson") {
    .haspar.poisson.Prior(obj, model.obj, verbose)
  } else if (dist == "binomial") {
    .haspar.binomial.Prior(obj, model.obj, verbose)
  } else if (dist == "exponential") {
    .haspar.exponential.Prior(obj, model.obj, verbose)
  } else if (dist == "cond.poisson") {
    .haspar.condpoisson.Prior(obj, model.obj, verbose)
  } else if (dist %in% c("normal", "normult")) {
    .haspar.normal.Prior(obj, model.obj, verbose)
  } else if (dist %in% c("student", "studmult")) {
    .haspar.student.Prior(obj, model.obj, verbose)
  }
}

### hasPar Prior Poisson
#' Checks for parameters in `prior` object for Poisson mixture
#' 
#' @description 
#' For internal usage only. This function checks, if a given 
#' [prior][prior-class] contains specified parameters in its slot `@@par`. For 
#' a Poisson mixture the parameters must be a `list` with named elements `a` 
#' and `b` for the Gamma shape and rate parameters, respectively. In addition 
#' the dimension of the parameters are checked for validity. The dimension of 
#' these parameters must conform to the number of components `K`. 
#' 
#' Hierarchical Poisson priors also need named shape parameter `g` and rate 
#' parameter `G` in the parameter list in slot `@@par`. The heorarchical prior 
#' is the same for each component, henceforth, there will be only one pair of 
#' such parameters.
#' 
#' @param obj A `prior` object to be checked.
#' @param model.obj A `model` object providing the model distribution for 
#'   which prior parameters should be checked. 
#' @param verbose A logical indicating, if the output should be verbose or 
#'   silent.
#' @returns Either a logical indicating, if the passed-in `prior` object 
#'   contains defined parameters or verbose output. Throws an error, if the 
#'   checks do not pass. 
#' @noRd
".haspar.poisson.Prior" <- function(obj, model.obj, verbose) {
  K <- model.obj@K
  if (length(obj@par) == 0) {
    if (verbose) {
      stop("Slot @par in 'prior' object is empty.",
        call. = FALSE
      )
    } else {
      return(FALSE)
    }
  } else {
    if (!("a" %in% names(obj@par))) {
      if (verbose) {
        stop(paste("Wrong specification of slot @par ",
          "in 'prior' object. Poisson models ",
          "need Gamma shape parameters named ",
          "'a'.",
          sep = ""
        ), call. = FALSE)
      } else {
        return(FALSE)
      }
    } else {
      if (dim(obj@par$a)[2] != K) {
        if (verbose) {
          stop(paste("Wrong specifcation of slot @par ",
            "in 'prior' object. Slot 'K' in ",
            "'model' object does not match ",
            "dimension of prior parameters.",
            sep = ""
          ), call. = FALSE)
        } else {
          return(FALSE)
        }
      } else {
        if (!("b" %in% names(obj@par))) {
          if (verbose) {
            stop(paste("Wrong specification of slot @par ",
              "in 'prior' object. Poisson models ",
              "need Gamma rate parameters named ",
              "'b'.",
              sep = ""
            ), call. = FALSE)
          } else {
            return(FALSE)
          }
        } else {
          if (dim(obj@par$b)[2] != K) {
            if (verbose) {
              stop(paste("Wrong specifcation of slot @par ",
                "in 'prior' object. Slot 'K' in ",
                "'model' object does not match ",
                "dimension of prior parameters.",
                sep = ""
              ), call. = FALSE)
            } else {
              return(FALSE)
            }
          } else {
            if (obj@hier) {
              if (!("g" %in% names(obj@par))) {
                if (verbose) {
                  stop(paste("Wrong specification of slot @par ",
                    "in 'prior' object if slot 'hier' ",
                    "is set to TRUE. Hierarchical Poisson models ",
                    "need Gamma shape hyperparameter named ",
                    "'g'.",
                    sep = ""
                  ), call. = FALSE)
                } else {
                  return(FALSE)
                }
              } else {
                if (!("G" %in% names(obj@par))) {
                  if (verbose) {
                    stop(paste("Wrong specification of slot @par ",
                      "in 'prior' object if slot 'hier' ",
                      "is set to TRUE. Hierarchical Poisson models ",
                      "need Gamma rate hyperparameter named ",
                      "'G'.",
                      sep = ""
                    ), call. = FALSE)
                  } else {
                    return(FALSE)
                  }
                } else {
                  return(TRUE)
                }
              }
            } else {
              return(TRUE)
            }
          }
        }
      }
    }
  }
}

#' Checks for parameters in `prior` object for conditional Poisson mixture
#' 
#' @description 
#' For internal usage only. This function checks, if a given 
#' [prior][prior-class] contains specified parameters in its slot `@@par`. For 
#' a conditional Poisson mixture the parameters must be a `list` with named 
#' elements `Q` and `N`, the component means and observations, respectively. 
#' Furthermore, parameters `a` and `b` are needed for each of the `K` 
#' components, defining the parameters of the uniform priors. As a final 
#' parameter conditional Poisson mixtures need the standard deviation `s` for 
#' the Metropolis-Hastings proposal. 
#' 
#' @param obj A `prior` object to be checked.
#' @param model.obj A `model` object providing the model distribution for 
#'   which prior parameters should be checked. 
#' @param verbose A logical indicating, if the output should be verbose or 
#'   silent.
#' @returns Either a logical indicating, if the passed-in `prior` object 
#'   contains defined parameters or verbose output. Throws an error, if the 
#'   checks do not pass. 
#' @noRd
".haspar.condpoisson.Prior" <- function(obj, model.obj, verbose) {
  K <- model.obj@K
  if (length(obj@par) == 0) {
    if (verbose) {
      stop("Slot @par in 'prior' object is empty.",
        call. = FALSE
      )
    } else {
      return(FALSE)
    }
  } else {
    if (!("Q" %in% names(obj@par))) {
      if (verbose) {
        stop(paste("Wrong specification of slot @par ",
          "in 'prior' object. Poisson models ",
          "need means named ",
          "'Q'.",
          sep = ""
        ), call. = FALSE)
      } else {
        return(FALSE)
      }
    } else {
      if (length(obj@par$Q) != K) {
        if (verbose) {
          stop(paste("Wrong specifcation of slot @par ",
            "in 'prior' object. Slot @K in ",
            "'model' object does not match ",
            "dimension of prior parameters.",
            sep = ""
          ), .call = FALSE)
        } else {
          return(FALSE)
        }
      } else {
        if (!("N" %in% names(obj@par))) {
          if (verbose) {
            stop(paste("Wrong specification of slot 'par' ",
              "in 'prior' object. Poisson models ",
              "need a number of observations per ",
              "component named ",
              "'N'.",
              sep = ""
            ), call. = FALSE)
          } else {
            return(FALSE)
          }
        } else {
          if (length(obj@par$N) != K) {
            if (verbose) {
              stop(paste("Wrong specifcation of slot 'par' ",
                "in 'prior' object. Slot @K in ",
                "'model' object does not match ",
                "dimension of prior parameters.",
                sep = ""
              ), call. = FALSE)
            } else {
              return(FALSE)
            }
          } else {
            if (!("a" %in% names(obj@par))) {
              if (verbose) {
                stop(paste(
                  "Wrong specification of slot @par ",
                  "in 'prior' object. Conditional ",
                  "Poisson models need a uniform ",
                  "distribution parameter 'a' ",
                  "defining the interval [a, b]."
                ),
                call. = FALSE
                )
              } else {
                return(FALSE)
              }
            } else {
              if (!("b" %in% names(obj@par))) {
                if (verbose) {
                  stop(paste(
                    "Wrong specification of slot @par ",
                    "in 'prior' object. Conditional ",
                    "Poisson models need a uniform ",
                    "distribution parameter 'b' ",
                    "defining the interval [a, b]."
                  ),
                  .call = FALSE
                  )
                } else {
                  return(FALSE)
                }
              } else {
                if (!("s" %in% names(obj@par))) {
                  if (verbose) {
                    stop(paste(
                      "Wrong specification of slot @par ",
                      "in 'prior' object. Conditional ",
                      "Poisson models need a parameter ",
                      "'s' defining the standard deviation ",
                      "of the Metropolis proposal."
                    ),
                    call. = FALSE
                    )
                  } else {
                    return(FALSE)
                  }
                } else {
                  return(TRUE)
                }
              }
            }
          }
        }
      }
    }
  }
}

#' Checks for parameters in `prior` object for Poisson mixture
#' 
#' @description 
#' For internal usage only. This function checks, if a given 
#' [prior][prior-class] contains specified parameters in its slot `@@par`. For 
#' a Binomial mixture the parameters must be a `list` with named elements `a` 
#' and `b` holding the shape and rate parameters of the Beta prior for each 
#' component. 
#' 
#' @param obj A `prior` object to be checked.
#' @param model.obj A `model` object providing the model distribution for 
#'   which prior parameters should be checked. 
#' @param verbose A logical indicating, if the output should be verbose or 
#'   silent.
#' @returns Either a logical indicating, if the passed-in `prior` object 
#'   contains defined parameters or verbose output. Throws an error, if the 
#'   checks do not pass. 
#' @noRd
".haspar.binomial.Prior" <- function(obj, model.obj, verbose) {
  K <- model.obj@K
  if (!length(obj@par)) {
    if (verbose) {
      stop("Slot @par in 'prior' object is empty.",
        call. = FALSE
      )
    } else {
      return(FALSE)
    }
  } else {
    if (!("a" %in% names(obj@par))) {
      if (verbose) {
        stop(paste("Wrong specification of slot @par ",
          "in 'prior' object. Binomial models ",
          "need Beta shape parameters named ",
          "'a'.",
          sep = ""
        ), call. = FALSE)
      } else {
        return(FALSE)
      }
    } else {
      if (dim(obj@par$a)[2] != K) {
        if (verbose) {
          stop(paste("Wrong specification of slot @par ",
            "in 'prior' object. Slot @K in ",
            "'model' object does not match ",
            "dimension of prior parameters.",
            sep = ""
          ), call. = FALSE)
        } else {
          return(FALSE)
        }
      } else {
        if (!("b" %in% names(obj@par))) {
          if (verbose) {
            stop(paste("Wrong specification of slot @par ",
              "in 'prior' object. Binomial models ",
              "need Beta rate parameters named ",
              "'b'.",
              sep = ""
            ), call. = FALSE)
          } else {
            return(FALSE)
          }
        } else {
          if (dim(obj@par$b)[2] != K) {
            if (verbose) {
              stop(paste("Wrong specification of slot @par ",
                "in 'prior' object. Slot @K in ",
                "'model' object does not match ",
                "dimension of prior parameters.",
                sep = ""
              ), call. = FALSE)
            } else {
              return(FALSE)
            }
          } else {
            return(TRUE)
          }
        }
      }
    }
  }
}

#' Checks for parameters in `prior` object for exponential mixture
#' 
#' @description 
#' For internal usage only. This function checks, if a given 
#' [prior][prior-class] contains specified parameters in its slot `@@par`. For 
#' a exponential mixture the parameters must be a `list` with named elements 
#' `a` and `b` defining the shape and rate parameters for the Gamma prior for 
#' each of the `K` components.
#' 
#' @param obj A `prior` object to be checked.
#' @param model.obj A `model` object providing the model distribution for 
#'   which prior parameters should be checked. 
#' @param verbose A logical indicating, if the output should be verbose or 
#'   silent.
#' @returns Either a logical indicating, if the passed-in `prior` object 
#'   contains defined parameters or verbose output. Throws an error, if the 
#'   checks do not pass. 
#' @noRd
".haspar.exponential.Prior" <- function(obj, model.obj, verbose) {
  K <- model.obj@K
  if (!length(obj@par)) {
    if (verbose) {
      stop("Slot @par in 'prior' object is empty.", call. = FALSE)
    } else {
      return(FALSE)
    }
  } else {
    if (!("a" %in% names(obj@par))) {
      if (verbose) {
        stop(paste("Wrong specification of slot @par ",
          "in 'prior' object. Priors for ",
          "exponential models need Gamma ",
          "shape parameters named 'a'.",
          sep = ""
        ), call. = FALSE)
      } else {
        return(FALSE)
      }
    } else {
      if (dim(obj@par$a)[2] != K) {
        if (verbose) {
          stop(paste("Wrong specification of slot @par ",
            "in 'prior' object. Slot @K in ",
            "'model' object does not match dimension ",
            "of prior parameters.",
            sep = ""
          ),
          call. = FALSE
          )
        } else {
          return(FALSE)
        }
      }
      if (!("b" %in% names(obj@par))) {
        if (verbose) {
          stop(paste("Wrong specification of slot @par ",
            "in 'prior' object. Priors for ",
            "exponential models need Gamma rate ",
            "parameters named 'b'.",
            sep = ""
          ),
          call. = FALSE
          )
        } else {
          return(FALSE)
        }
      } else {
        if (dim(obj@par$b)[2] != K) {
          if (verbose) {
            stop(paste("Wrong specification of slot @par ",
              "in 'prior' object. Slot @K in ",
              "'model' object does not match dimension ",
              "of prior parameters.",
              sep = ""
            ),
            call. = FALSE
            )
          } else {
            return(FALSE)
          }
        } else {
          return(TRUE)
        }
      }
    }
  }
}

#' Checks for parameters in `prior` object for a normal mixture
#' 
#' @description 
#' For internal usage only. This function checks, if a given 
#' \code{\link{prior-class}} contains specified parameters in its slot `@@par`. For 
#' a normal mixture the parameters must be a `list` with named elements `mu`, 
#' and `sigma` defining the prior parameters for the mean and standard 
#' deviations respectively. 
#' 
#' ## Conditionally conjugate prior
#' In case a conditional conjugate prior is chosen 
#' `mu` and `sigma` must be lists with elements `b` and `N0` and `c` and `C`, 
#' respectively. `b` and `N0` define the parameters of a normal prior with 
#' means `b` and standard deviations `sigma/N0`. `c` and `C` define the shape 
#' and rate parameters of an inverse Gamma prior for the standard deviations.
#' 
#' ## Independent prior
#' If an independent prior is chosen, the elements `mu` and `sigma` in slot 
#' `@@par` of the `prior` object must contain lists with the following 
#' elements: `b` and `Binv` and `c` and `C`. `b` and `Binv` are the means and 
#' inverse standard deviations of normal priors and `c` and `C` are the shape 
#' and rate parameters of an inverse Gamma distribution, respectively. 
#' 
#' ## Hierarchical prior
#' In case of an hierarchical prior the list referred to by the name `sigma` in 
#' `@@par` needs to contain a shape parameter `g` and rate parameter `G` of the 
#' hierarchical gamma prior for the prior parameter `C`. 
#' 
#' ## Multivariate normal mixtures
#' In case that the `model.obj` defines a multivariate normal mixture 
#' distribution the prior parameters are defined by `list` with elements `mu` 
#' and `sigma` for the prior for the means and the prior for the 
#' variance-covariance matrices, respectively. 
#' 
#' ### Conditionally conjugate prior 
#' In case of a conditionally conjugate prior, the prior for the means is 
#' defined by a `list` with elements named `b` and `N0` with `b` an `rxK` 
#' matrix defining the means of the normal prior and `N0` an `1xK` vector 
#' defining the scaling constants for the standard deviations of the normal 
#' prior for the means. 
#' The element `sigma` is a `list` with elements `c` and `C` defining the 
#' parameters of the Wishart prior for the covariance matrices. `c` is an 
#' `1xK` matrix or vector and `C` an `rxrxK` array. In 
#' addition an element `logdetC` is required that contains the logarithmized 
#' determinants of the matrices in `C`. 
#' 
#' ### Independent prior
#' In case an independent prior is used the element `mu` must contain a `list` 
#' with elements `b` and `Binv`. `b` is the means matrix of dimension 
#' `rxK` for the normal prior of `mu` and `Binv` contain the inverted 
#' variance-covariance matrices of the normal prior. 
#' The corresponding prior for the variance-covariance matrices `sigma` is the 
#' same as for the conditionally conjugate prior.
#' 
#' @param obj A `prior` object to be checked.
#' @param model.obj A `model` object providing the model distribution for 
#'   which prior parameters should be checked. 
#' @param verbose A logical indicating, if the output should be verbose or 
#'   silent.
#' @returns Either a logical indicating, if the passed-in `prior` object 
#'   contains defined parameters or verbose output. Throws an error, if the 
#'   checks do not pass. 
#' @noRd
".haspar.normal.Prior" <- function(obj, model.obj, verbose) {
  K <- model.obj@K
  if (length(obj@par) == 0) {
    if (verbose) {
      stop("Slot @par in 'prior' object is empty.",
        call. = FALSE
      )
    } else {
      return(FALSE)
    }
  } else {
    if (!"mu" %in% names(obj@par)) {
      if (verbose) {
        stop(paste("Wrong specification of slot @par in ",
          "'prior' object. Priors for Normal models ",
          "need a 'list' object called 'mu'.",
          sep = ""
        ), call. = FALSE)
      } else {
        return(FALSE)
      }
    } else {
      if (!"b" %in% names(obj@par$mu)) {
        if (verbose) {
          stop(paste("Wrong specification of slot @par in ",
            "'prior' object. Priors for Normal models ",
            "need mean values named 'b'.",
            sep = ""
          ),
          call. = FALSE
          )
        } else {
          return(FALSE)
        }
      } else {
        if (obj@type == "condconjugate") {
          if (!"N0" %in% names(obj@par$mu)) {
            if (verbose) {
              stop(paste("Wrong specification of slot @par in ",
                "'prior' object. Conditionally conjugate ",
                "priors for Normal models need variances ",
                "named 'N0'.",
                sep = ""
              ), call. = FALSE)
            } else {
              return(FALSE)
            }
          } else {
            if (!"sigma" %in% names(obj@par)) {
              if (verbose) {
                stop(paste("Wrong specification of slot @par in ",
                  "'prior' object. Priors for Normal models ",
                  "need a 'list' object named 'sigma'.",
                  sep = ""
                ),
                call. = FALSE
                )
              } else {
                return(FALSE)
              }
            } else {
              if (!"c" %in% names(obj@par$sigma)) {
                if (verbose) {
                  stop(paste("Wrong specification of slot @par in ",
                    "'prior' object. Priors for Normal models ",
                    "need a inverse Gamma shape parameter named ",
                    "'c'.",
                    sep = ""
                  ), call. = FALSE)
                } else {
                  return(FALSE)
                }
              } else {
                if (!"C" %in% names(obj@par$sigma)) {
                  if (verbose) {
                    stop(paste("Wrong specification of slot @par in ",
                      "'prior' object'. Priors for Normal models ",
                      "need a inverse Gamma rate parameter named ",
                      "'C'.",
                      sep = ""
                    ), call. = FALSE)
                  }
                } else {
                  if (model.obj@r > 1 && !"logdetC" %in% names(obj@par$sigma)) {
                    if (verbose) {
                      stop(paste("Wrong specification of slot @par in ",
                        "'prior' object. Priors for multivariate ",
                        "Normal models need the logarithmised ",
                        "determinant of the shape matrix of the ",
                        "Wishart distribution and has to be named ",
                        "'logdetC'.",
                        sep = ""
                      ), call. = FALSE)
                    } else {
                      return(FALSE)
                    }
                  } else {
                    if (obj@hier && !"g" %in% names(obj@par$sigma)) {
                      if (verbose) {
                        stop(paste("Wrong specification of slot @par in ",
                          "'prior' object. Hierarchical priors for ",
                          "Normal models need a Gamma shape parameter ",
                          "named 'g'.",
                          sep = ""
                        ), call. = FALSE)
                      } else {
                        return(FALSE)
                      }
                    } else if (obj@hier && !"G" %in% names(obj@par$sigma)) {
                      if (verbose) {
                        stop(paste("Wrong specification of slot @par in ",
                          "'prior' object. Hierarchical priors for ",
                          "Normal models need a Gamma rate parameter ",
                          "named 'G'.",
                          sep = ""
                        ), call. = FALSE)
                      } else {
                        return(FALSE)
                      }
                    } else {
                      return(TRUE)
                    }
                  }
                }
              }
            }
          }
        } else {
          if (!"Binv" %in% names(obj@par$mu)) {
            if (verbose) {
              stop(paste("Wrong specification of slot @par in ",
                "'prior' object. Independent priors for ",
                "Normal models need variances named 'Binv'.",
                sep = ""
              ), call. = FALSE)
            } else {
              return(FALSE)
            }
          } else {
            if (!"sigma" %in% names(obj@par)) {
              if (verbose) {
                stop(paste("Wrong specification of slot @par in ",
                  "'prior' object. Priors for Normal models ",
                  "need a 'list' object named 'sigma'.",
                  sep = ""
                ),
                call. = FALSE
                )
              } else {
                return(FALSE)
              }
            } else {
              if (!"c" %in% names(obj@par$sigma)) {
                if (verbose) {
                  stop(paste("Wrong specification of slot @par in ",
                    "'prior' object. Priors for Normal models ",
                    "need a inverse Gamma shape parameter named ",
                    "'c'.",
                    sep = ""
                  ), call. = FALSE)
                } else {
                  return(FALSE)
                }
              } else {
                if (!"C" %in% names(obj@par$sigma)) {
                  if (verbose) {
                    stop(paste("Wrong specification of slot @par in ",
                      "'prior' object. Priors for Normal models ",
                      "need a inverse Gamma rate parameter named ",
                      "'C'.",
                      sep = ""
                    ), call. = FALSE)
                  } else {
                    return(FALSE)
                  }
                } else {
                  if (model.obj@r > 1 && !"logdetC" %in% names(obj@par$sigma)) {
                    if (verbose) {
                      stop(paste("Wrong specification of slot @par in ",
                        "'prior' object. Priors for multivariate ",
                        "Normal models need the logarithmised ",
                        "determinant of the shape matrix of the ",
                        "Wishart distribution and has to be named ",
                        "'logdetC'.",
                        sep = ""
                      ), call. = FALSE)
                    } else {
                      return(FALSE)
                    }
                  } else {
                    if (obj@hier && !"g" %in% names(obj@par$sigma)) {
                      if (verbose) {
                        stop(paste("Wrong specification of slot @par in ",
                          "'prior' object. Hierarchical priors for ",
                          "Normal models need a Gamma shape parameter ",
                          "named 'g'.",
                          sep = ""
                        ), call. = FALSE)
                      } else {
                        return(FALSE)
                      }
                    } else if (obj@hier && !"G" %in% names(obj@par$sigma)) {
                      if (verbose) {
                        stop(paste("Wrong specification of slot @par in ",
                          "'prior' object. Hierarchical priors for ",
                          "Normal models need a Gamma rate parameter ",
                          "named 'G'.",
                          sep = ""
                        ), call. = FALSE)
                      } else {
                        return(FALSE)
                      }
                    } else {
                      return(TRUE)
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

#' Checks for parameters in `prior` object for a normal mixture
#' 
#' @description 
#' For internal usage only. This function checks, if a given 
#' [prior][prior-class] contains specified parameters in its slot `@@par`. For 
#' a normal mixture the parameters must be a `list` with named elements `mu`, 
#' and `sigma` defining the prior parameters for the mean and standard 
#' deviations respecitvely. 
#' 
#' ## Conditionally conjugate prior
#' In case a conditional conjugate prior is chosen 
#' `mu` and `sigma` must be lists with elements `b` and `N0` and `c` and `C`, 
#' respectively. `b` and `N0` define the parameters of a normal prior with 
#' means `b` and standard deviations `sigma/N0`. `c` and `C` define the shape 
#' and rate parameters of an inverse Gamma prior for the standard deviations.
#' 
#' ## Independent prior
#' If an independent prior is chosen, the elements `mu` and `sigma` in slot 
#' `@@par` of the `prior` object must contain lists with the following 
#' elements: `b` and `Binv` and `c` and `C`. `b` and `Binv` are the means and 
#' inverse standard deviations of normal priors and `c` and `C` are the shape 
#' and rate parameters of an inverse Gamma distribution respectively. 
#' 
#' ## Hierarchical prior
#' In case of an hierarchical prior the list referred to by the name `sigma` in 
#' `@@par` needs to contain a shape parameter `g` and rate parameter `G` of the 
#' hierarchical gamma prior for the prior parameter `C`. 
#' 
#' ## Multivariate normal mixtures
#' In case that the `model.obj` defines a multivariate normal mixture 
#' distribution the prior parameters are defined by `list` with elements `mu` 
#' and `sigma` for the prior for the means and the prior for the 
#' variance-covariance matrices, resepectively. 
#' 
#' ### Conditionally conjugate prior 
#' In case of a conditionally conjugate prior, the prior for the means is 
#' defined by a `list` with elements named `b` and `N0` with `b` an `rxK` 
#' matrix defining the means of the normal prior and `N0` an `1xK` vector 
#' defining the scaling constants for the standard deviations of the normal 
#' prior for the means. 
#' The element `sigma` is a `list` with elements `c` and `C` defining the 
#' parameters of the Wishart prior for the covariance matrices. `c` is an 
#' `1xK` matrix or vector and `C` an `rxrxK` array. In 
#' addition an element `logdetC` is required that contains the logarithmized 
#' determinants of the matrices in `C`. 
#' 
#' ### Independent prior
#' In case an independent prior is used the element `mu` must contain a `list` 
#' with elements `b` and `Binv`. `b` is the means matrix of dimension 
#' `rxK` for the normal prior of `mu` and `Binv` contain the inverted 
#' variance-covariance matrices of the normal prior. 
#' The corresponding prior for the variance-covariance matrices `sigma` is the 
#' same as for the conditionally conjugate prior.
#' 
#' ## Prior for the degrees of freedom
#' The prior for the degrees of freedom is the same for univariate and 
#' multivariate mixtures. in both cases the slot `@@par` must contain in its 
#' `list` of prior parameters a field named `df` that further contains a `list` 
#' with elements `type`, `trans`, `a0`, `b0`, and `d`. `type` defines the type 
#' of prior used and must be set to `"inhier"` for the independent hierarchcial 
#' prior. For this prior `trans` is the translation parameter and the other 
#' parameters are further parameters to define the prior. All parameters are 
#' `numeric`s. Furthermore, there is an additional parameter named `mhtune` 
#' defining the width parameters of the uniform log random walk proposals for 
#' the degrees of freedom in Metropolis-Hastings sampling. This field has to be 
#' a vector of size `1xK`.
#' 
#' @param obj A `prior` object to be checked.
#' @param model.obj A `model` object providing the model distribution for 
#'   which prior parameters should be checked. 
#' @param verbose A logical indicating, if the output should be verbose or 
#'   silent.
#' @returns Either a logical indicating, if the passed-in `prior` object 
#'   contains defined parameters or verbose output. Throws an error, if the 
#'   checks do not pass. 
#' @noRd
".haspar.student.Prior" <- function(obj, model.obj, verbose) {
  K <- model.obj@K
  if (length(obj@par) == 0) {
    if (verbose) {
      stop("Slot @par in 'prior' object is empty.",
        call. = FALSE
      )
    } else {
      return(FALSE)
    }
  } else {
    if (!"mu" %in% names(obj@par)) {
      if (verbose) {
        stop(paste("Wrong specification of slot @par in ",
          "'prior' object. Priors for Student-t models ",
          "need a 'list' object called 'mu'.",
          sep = ""
        ), call. = FALSE)
      } else {
        return(FALSE)
      }
    } else {
      if (!"b" %in% names(obj@par$mu)) {
        if (verbose) {
          stop(paste("Wrong specification of slot @par in ",
            "'prior' object. Priors for Student-t models ",
            "need mean values named 'b'.",
            sep = ""
          ),
          call. = FALSE
          )
        } else {
          return(FALSE)
        }
      } else {
        if (obj@type == "condconjugate") {
          if (!"N0" %in% names(obj@par$mu)) {
            if (verbose) {
              stop(paste("Wrong specification of slot @par in ",
                "'prior' object. Conditionally conjugate ",
                "priors for Student-t models need variances ",
                "named 'N0'.",
                sep = ""
              ), call. = FALSE)
            } else {
              return(FALSE)
            }
          } else {
            if (!"sigma" %in% names(obj@par)) {
              if (verbose) {
                stop(paste("Wrong specification of slot @par in ",
                  "'prior' object. Priors for Student-t models ",
                  "need a 'list' object named 'sigma'.",
                  sep = ""
                ),
                call. = FALSE
                )
              } else {
                return(FALSE)
              }
            } else {
              if (!"c" %in% names(obj@par$sigma)) {
                if (verbose) {
                  stop(paste("Wrong specification of slot @par in ",
                    "'prior' object. Priors for Student-t models ",
                    "need a inverse Gamma shape parameter named ",
                    "'c'.",
                    sep = ""
                  ), call. = FALSE)
                } else {
                  return(FALSE)
                }
              } else {
                if (!"C" %in% names(obj@par$sigma)) {
                  if (verbose) {
                    stop(paste("Wrong specification of slot @par in ",
                      "'prior' object'. Priors for Student-t models ",
                      "need a inverse Gamma rate parameter named ",
                      "'C'.",
                      sep = ""
                    ), call. = FALSE)
                  }
                } else {
                  if (model.obj@r > 1 && !"logdetC" %in% names(obj@par$sigma)) {
                    if (verbose) {
                      stop(paste("Wrong specification of slot @par in ",
                        "'prior' object. Priors for multivariate ",
                        "Student-t models need the logarithmised ",
                        "determinant of the shape matrix of the ",
                        "Wishart distribution and has to be named ",
                        "'logdetC'.",
                        sep = ""
                      ), call. = FALSE)
                    } else {
                      return(FALSE)
                    }
                  } else {
                    if (obj@hier && !"g" %in% names(obj@par$sigma)) {
                      if (verbose) {
                        stop(paste("Wrong specification of slot @par in ",
                          "'prior' object. Hierarchical priors for ",
                          "Student-t models need a Gamma shape parameter ",
                          "named 'g'.",
                          sep = ""
                        ), call. = FALSE)
                      } else {
                        return(FALSE)
                      }
                    } else if (obj@hier && !"G" %in% names(obj@par$sigma)) {
                      if (verbose) {
                        stop(paste("Wrong specification of slot @par in ",
                          "'prior' object. Hierarchical priors for ",
                          "Student-t models need a Gamma rate parameter ",
                          "named 'G'.",
                          sep = ""
                        ), call. = FALSE)
                      } else {
                        return(FALSE)
                      }
                    } else {
                      if (!"df" %in% names(obj@par)) {
                        if (verbose) {
                          stop(paste("Wrong specification of slot @par in ",
                            "'prior' object. Priors for Student-t ",
                            "models need a 'list' object named ",
                            "'df'.",
                            sep = ""
                          ), call. = FALSE)
                        } else {
                          return(FALSE)
                        }
                      } else {
                        if (!"type" %in% names(obj@par$df)) {
                          if (verbose) {
                            stop(paste("Wrong specification of slot @par in ",
                              "'prior' object. Priors for the degrees of ",
                              "freedom in Student-t models need a type ",
                              "named 'type'.",
                              sep = ""
                            ), call. = FALSE)
                          } else {
                            return(FALSE)
                          }
                        } else {
                          if (!all(c("trans", "a0", "b0", "d") %in% names(obj@par$df))
                          ) {
                            if (verbose) {
                              stop(paste("Wrong specification of slot @par in ",
                                "'prior' object. Priors for the degrees of ",
                                "freedom in Student-t models need ",
                                "hyperparameters named 'trans', 'a0', 'b0' ",
                                "and 'd'.",
                                sep = ""
                              ), call. = FALSE)
                            } else {
                              return(FALSE)
                            }
                          } else {
                            if (!"mhtune" %in% names(obj@par$df)) {
                              if (verbose) {
                                stop(paste("Wrog specification of slot @par in ",
                                  "'prior' object. Priors for the degrees ",
                                  "of freedom need Metropolis-Hastings ",
                                  "tuning parameters named 'mhtune'.",
                                  sep = ""
                                ), call. = FALSE)
                              } else {
                                return(FALSE)
                              }
                            } else {
                              return(TRUE)
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        } else {
          if (!"Binv" %in% names(obj@par$mu)) {
            if (verbose) {
              stop(paste("Wrong specification of slot @par in ",
                "'prior' object. Independent priors for ",
                "Normal models need variances named 'Binv'.",
                sep = ""
              ), call. = FALSE)
            } else {
              return(FALSE)
            }
          } else {
            if (!"sigma" %in% names(obj@par)) {
              if (verbose) {
                stop(paste("Wrong specification of slot @par in ",
                  "'prior' object. Priors for Normal models ",
                  "need a 'list' object named 'sigma'.",
                  sep = ""
                ),
                call. = FALSE
                )
              } else {
                return(FALSE)
              }
            } else {
              if (!"c" %in% names(obj@par$sigma)) {
                if (verbose) {
                  stop(paste("Wrong specification of slot @par in ",
                    "'prior' object. Priors for Normal models ",
                    "need a inverse Gamma shape parameter named ",
                    "'c'.",
                    sep = ""
                  ), call. = FALSE)
                } else {
                  return(FALSE)
                }
              } else {
                if (!"C" %in% names(obj@par$sigma)) {
                  if (verbose) {
                    stop(paste("Wrong specification of slot @par in ",
                      "'prior' object. Priors for Normal models ",
                      "need a inverse Gamma rate parameter named ",
                      "'C'.",
                      sep = ""
                    ), call. = FALSE)
                  } else {
                    return(FALSE)
                  }
                } else {
                  if (model.obj@r > 1 && !"logdetC" %in% names(obj@par$sigma)) {
                    if (verbose) {
                      stop(paste("Wrong specification of slot @par in ",
                        "'prior' object. Priors for multivariate ",
                        "Normal models need the logarithmised ",
                        "determinant of the shape matrix of the ",
                        "Wishart distribution and has to be named ",
                        "'logdetC'.",
                        sep = ""
                      ), call. = FALSE)
                    } else {
                      return(FALSE)
                    }
                  } else {
                    if (obj@hier && !"g" %in% names(obj@par$sigma)) {
                      if (verbose) {
                        stop(paste("Wrong specification of slot @par in ",
                          "'prior' object. Hierarchical priors for ",
                          "Normal models need a Gamma shape parameter ",
                          "named 'g'.",
                          sep = ""
                        ), call. = FALSE)
                      } else {
                        return(FALSE)
                      }
                    } else if (obj@hier && !"G" %in% names(obj@par$sigma)) {
                      if (verbose) {
                        stop(paste("Wrong specification of slot @par in ",
                          "'prior' object. Hierarchical priors for ",
                          "Normal models need a Gamma rate parameter ",
                          "named 'G'.",
                          sep = ""
                        ), call. = FALSE)
                      } else {
                        return(FALSE)
                      }
                    } else {
                      if (!"df" %in% names(obj@par)) {
                        if (verbose) {
                          stop(paste("Wrong specification of slot @par in ",
                            "'prior' object. Priors for Student-t ",
                            "models need a 'list' object named ",
                            "'df'.",
                            sep = ""
                          ), call. = FALSE)
                        } else {
                          return(FALSE)
                        }
                      } else {
                        if (!"type" %in% names(obj@par$df)) {
                          if (verbose) {
                            stop(paste("Wrong specification of slot @par in ",
                              "'prior' object. Priors for the degrees of ",
                              "freedom in Student-t models need a type ",
                              "named 'type'.",
                              sep = ""
                            ), call. = FALSE)
                          } else {
                            return(FALSE)
                          }
                        } else {
                          if (!all(c("trans", "a0", "b0", "d") %in% names(obj@par$df))
                          ) {
                            if (verbose) {
                              stop(paste("Wrong specification of slot @par in ",
                                "'prior' object. Priors for the degrees of ",
                                "freedom in Student-t models need ",
                                "hyperparameters named 'trans', 'a0', 'b0' ",
                                "and 'd'.",
                                sep = ""
                              ), call. = FALSE)
                            } else {
                              return(FALSE)
                            }
                          } else {
                            if (!"mhtune" %in% names(obj@par$df)) {
                              if (verbose) {
                                stop(paste("Wrong specification of slot @par in ",
                                  "'prior' object. Priors for the degrees ",
                                  "of freedom need Metropolis-Hastings ",
                                  "tuning parameters named 'mhtune'.",
                                  sep = ""
                                ), call. = FALSE)
                              } else {
                                return(FALSE)
                              }
                            } else {
                              return(TRUE)
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

#' Generate default prior for Poisson mixture
#' 
#' @description 
#' For internal usage only. This function constructs the default priors for 
#' a Poisson mixture. 
#' 
#' @details 
#' The type of a data-dependent Poisson prior is
#'                 is always conditionally conjugate Gamma with
#'                 parameters:     a: shape,   1 x model.obj@K
#'                 b: rate,    1 x model.obj@K
#'                 If not otherwise specified in 'varargin' an
#'                 hierarchical Gamma distribution is chosen with
#'                 parameters:     g: shape,   1 x 1
#'                                 G: rate,    1 x 1.
#' 
#' @param obj A `prior` object to be enriched by prior parameters. 
#' @param fdata.obj An `fdata` object containing the data.
#' @param model.obj A `model` object specifying the mixture model.
#' @param varargin A `prior` object passed in by the user. Optional.
#' @return A `prior` object with specified prior parameters for the prior of 
#'   a Poisson mixture model.
#' @noRd
#' @seealso 
#' * [priordefine()] for the calling function.
".generatePriorPoisson" <- function(obj, fdata.obj, model.obj,
                                    varargin) {
  K <- model.obj@K
  datam <- getColY(fdata.obj)
  if (is.null(varargin)) {
    obj@hier <- TRUE ## default prior is hierarchical
    obj@type <- "condconjugate"
  } else {
    obj@hier <- varargin@hier
    obj@type <- "condconjugate"
  }
  ## Default prior based on matching moments
  ## See Frühwirth-Schnatter (2006) and Viallefont et. al. (2002)
  ## for more information.

  ## Choose level of overdispersion, depending on the
  ## ratio overdispersion/mean^2
  ## no idea: data-based choice
  mean <- mean(datam, na.rm = TRUE)
  over <- as.numeric(var(datam, na.rm = TRUE) - mean)
  if (over > 0) {
    a0 <- mean^2 / over
  } else {
    a0 <- 10
  }
  if (obj@hier) {
    g0 <- 0.5
    G0 <- mean * g0 / a0
    b0 <- g0 / G0
    par <- list(
      a = array(a0, dim = c(1, K)),
      b = array(b0, dim = c(1, K)),
      g = g0, G = G0
    )
  } else {
    b0 <- a0 / mean
    par <- list(
      a = array(a0, dim = c(1, K)),
      b = array(b0, dim = c(1, K))
    )
  }
  obj@par <- par
  return(obj)
}

#' Generate default prior for conditional Poisson mixture
#' 
#' @description 
#' For internal usage only. This function constructs the default priors for 
#' a conditional Poisson mixture. 
#' 
#' @param obj A `prior` object to be enriched by prior parameters. 
#' @param fdata.obj An `fdata` object containing the data.
#' @param model.obj A `model` object specifying the mixture model.
#' @param varargin A `prior` object passed in by the user. Optional.
#' @return A `prior` object with specified prior parameters for the prior of 
#'   a conditional Poisson mixture model.
#' @importFrom stats kmeans
#' @noRd
#' @seealso 
#' * [priordefine()] for the calling function.
".generatePriorCondPoisson" <- function(obj, fdata.obj, model.obj,
                                        s) {
  K <- model.obj@K
  if (fdata.obj@bycolumn) {
    datam <- fdata.obj@y
  } else {
    datam <- t(fdata.obj@y)
  }
  obj@hier <- FALSE ## default prior is hierarchical
  obj@type <- "condconjugate"
  clust <- kmeans(datam, model.obj@K)
  Q <- vector("numeric", K)
  N <- vector("numeric", K)
  for (k in 1:K) {
    Q[k] <- mean(fdata.obj@y[clust$cluster == k], na.rm = TRUE)
    N[k] <- sum(clust$cluster == k)
  }
  Qsort <- sort(Q, decreasing = FALSE, index.return = TRUE)
  Q <- Qsort$x
  N <- N[Qsort$ix]
  ## Q   <- mean( datam, na.rm = TRUE )
  ## if ( nrow( datam ) %% K == 0 ) {
  ##    N   <- nrow( datam ) / K
  ## } else {
  ##    N <- nrow( datam ) / K
  ##    N <- c( ceiling( N ), floor( N ) )
  ## }
  mu1 <- Q[1]
  mu2 <- mean(datam[clust$cluster == which(Qsort$ix == 1)]^2, na.rm = TRUE)
  A <- max(0, mu1 - sqrt(3 * (mu2 - mu1^2)))
  B <- mu1 + sqrt(3 * (mu2 - mu1^2))
  pars <- list(
    Q = as.array(Q),
    N = as.array(N),
    a = A,
    b = B,
    s = s
  )
  obj@par <- pars
  return(obj)
}

#' Selects the Beta prior by quantiles. 
#' 
#' @description 
#' For internal usage only. Not used. Rather a relic. 
#'
#' @param quantile1 A `list`. 
#' @param quantile2 A `list`.
#' @return Unknown.
#' @importFrom stats pbeta approx
#' @noRd
".select.beta.Prior" <- function(quantile1, quantile2) {
  betaprior1 <- function(K, x, p) {
    m.lo <- 0.0
    m.hi <- 1
    flag <- 0
    while (flag == 0) {
      m0 <- (m.lo + m.hi) / 2
      p0 <- pbeta(x, K * m0, K * (1 - m0))
      if (p0 < p) m.hi <- m0 else m.lo <- m0
      if (abs(p0 - p) < .0001) flag <- 1
    }
    return(m0)
  }
  p1 <- quantile1$p
  x1 <- quantile1$x
  p2 <- quantile2$p
  x2 <- quantile2$x
  logK <- seq(-3, 8, length = 100)
  K <- exp(logK)
  m <- sapply(K, betaprior1, x1, p1)
  prob2 <- pbeta(x2, K * m, K * (1 - m))
  ind <- ((prob2 > 0) & (prob2 < 1))
  app <- approx(prob2[ind], logK[ind], p2)
  K0 <- exp(app$y)
  m0 <- betaprior1(K0, x1, p1)
  return(round(K0 * c(m0, (1 - m0)), 2))
}

#' Generate default prior for Binomial mixture
#' 
#' @description 
#' For internal usage only. This function constructs the default priors for 
#' a Binomial mixture. 
#' 
#' @details 
#' The type of generated Binomial prior is always a conditionally conjugate, 
#' i.e. Beta with parameters:                
#' a:  shape,  1 x model.obj@@K
#' b:  rate,   1 x model.obj@@K;
#' starting values are a = (1, 1), b = (1, 1).
#' 
#' @param obj A `prior` object to be enriched by prior parameters. 
#' @param model.obj A `model` object specifying the mixture model.
#' @return A `prior` object with specified prior parameters for the prior of 
#'   a Binomial mixture model.
#' @noRd
#' @seealso 
#' * [priordefine()] for the calling function.
".generatePriorBinomial" <- function(obj, model.obj) {
  K <- model.obj@K
  obj@type <- "condconjugate"
  ## uniform prior ##
  a0 <- 1
  b0 <- 1
  obj@par <- list(
    a = array(a0, dim = c(1, K)),
    b = array(b0, dim = c(1, K))
  )
  return(obj)
}

#' Generate default prior for exponential mixture
#' 
#' @description 
#' For internal usage only. This function constructs the default priors for 
#' a exponential mixture. 
#' 
#' @details 
#' If the identifier `prior.wagner == TRUE` the prior from Wagner (2007) is 
#' taken. In the remaining case a prior is constructed from the analysis of 
#' over-dispersion in the observations. This prior can also be hierarchical
#' if specified.
#' 
#' @param obj A `prior` object to be enriched by prior parameters. 
#' @param model.obj A `model` object specifying the mixture model.
#' @return A `prior` object with specified prior parameters for the prior of 
#'   a exponential mixture model.
#' @noRd
#' @seealso 
#' * [priordefine()] for the calling function.
".generatePriorExponential" <- function(obj, fdata.obj, model.obj,
                                        varargin, prior.wagner) {
  if (is.null(varargin)) {
    obj@hier <- TRUE
  } else {
    obj@hier <- varargin@hier
  }
  obj@type <- "condconjugate"
  datam <- getColY(fdata.obj)
  K <- model.obj@K
  if (prior.wagner) {
    # Prior following Wagner (2007) ##
    obj@hier <- FALSE
    a0 <- 0.1
    be <- mean(datam, na.rm = TRUE) * a0
    obj@par <- list(
      a = array(a0, dim = c(1, K)),
      b = array(be, dim = c(1, K))
    )
  } else {
    # Prior based on matching moments
    # Choose level of overdispersion
    # levover
    #   0:  no idea, data based choice
    #   1:  low degree of overdispersion; ratio smaller than 1
    #   2:  medium degree of overdispersion; ratio close to 1
    #   3:  high degree of overdispersion; ratio larger than 1
    levover <- 0
    a00 <- c(10, 8 / 3, 2.1)
    d.mean <- mean(datam, na.rm = TRUE)
    d2.mean <- mean(datam^2, na.rm = TRUE)
    over <- sd(datam, na.rm = TRUE) / d.mean - 1
    if (levover == 0) {
      if (d2.mean - 2 * d.mean^2 > 0) {
        a0 <- 2 * var(datam) / (d2.mean - 2 * d.mean^2)
      } else {
        if (over < 0.2) {
          levover <- 1
          a0 <- a00[1]
        } else if (over < 2) {
          levover <- 2
          a0 <- a00[2]
        } else {
          levover <- 3
          a0 <- a00[3]
        }
      }
    } else {
      a0 <- a00(levover)
    }
    if (obj@hier) {
      g0 <- 0.5
      G0 <- g0 / d.mean * (a0 - 1)
      be <- g0 / G0
      obj@par <- list(
        a = array(a0, dim = c(1, K)),
        b = array(be, dim = c(1, K)),
        g = g0,
        G = G0
      )
    } else {
      be <- d.mean * (a0 - 1)
      obj@par <- list(
        a = array(a0, dim = c(1, K)),
        b = array(be, dim = c(1, K))
      )
    }
  }
  return(obj)
}

#' Generate default prior for normal or Student-t mixtures
#' 
#' @description 
#' For internal usage only. This function constructs the default priors for 
#' a normal or Student-t mixture. 
#' 
#' @param obj A `prior` object to be enriched by prior parameters. 
#' @param data.obj An `fdata` object storing the observations.
#' @param model.obj A `model` object specifying the mixture model.
#' @param varargin A `rpior` object passed in by the user with predefined 
#'   slots.
#' @return A `prior` object with specified prior parameters for the prior of 
#'   a normal or Student-t mixture model.
#' @noRd
#' @seealso 
#' * [priordefine()] for the calling function.
".generatePriorNorstud" <- function(obj, data.obj,
                                    model.obj, varargin) {
  r <- data.obj@r
  K <- model.obj@K
  datam <- getColY(data.obj)
  ## check if varargin is non-empty and prior object ##
  ## set hierarchical or non-hierarchical prior ##
  if (is.null(varargin)) {
    ## default prior: independent hierarchical prior ##
    obj@hier <- hier <- TRUE
    obj@type <- "independent"
  } else {
    obj@hier <- hier <- varargin@hier
    obj@type <- varargin@type
  }
  conjugate.prior <- obj@type == "condconjugate"
  bensmail <- FALSE
  rich.green <- FALSE
  if (conjugate.prior || !hier) {
    bensmail <- TRUE ## Prior following Bensmail et al.
  } else {
    rich.green <- TRUE ## Prior following Richardson and Green for r = 1
    ##                 Stephens (1997a)     for r = 2 only
  }
  if (rich.green) { ## Richardson and Green (1997) or Stephens (1997a)
    ## row vectors: dimension 1 x r
    max <- apply(datam, 2, max, na.rm = TRUE)
    min <- apply(datam, 2, min, na.rm = TRUE)
    mean <- (max + min) * .5
    cov <- diag((max - min)^2, nrow = model.obj@r)
  } else {
    ## row vectors: dimension 1 x r
    mean <- apply(datam, 2, mean, na.rm = TRUE)
    cov <- var(datam, na.rm = TRUE)
  }
  b0 <- mean
  if (conjugate.prior) {
    B0sc <- 1 ## info contained in a standard conjugate (sc) prior (equal to N0)
    ## Bensmail et al. (1997)
  } else {
    B0inv <- solve(cov) ## info contained in a non-conjugate prior,
    ## i.e. either by Richardson Green (1997)
  }
  if (!conjugate.prior) {
    if (r > 1) {
      par.mu <- list(
        b = array(t(b0), dim = c(r, K)),
        Binv = array(B0inv, dim = c(r, r, K))
      )
    } else { ## r = 1
      par.mu <- list(
        b = array(b0, dim = c(1, K)),
        Binv = array(B0inv, dim = c(1, K))
      )
    }
  } else { ## conditionally conjugate prior
    if (r > 1) {
      par.mu <- list(
        b = array(t(b0), dim = c(r, K)),
        N0 = array(B0sc, dim = c(1, K))
      )
    } else { ## r = 1
      par.mu <- list(
        b = array(b0, dim = c(1, K)),
        N0 = array(B0sc, dim = c(1, K))
      )
    }
  }

  ## prior sigma ##
  ## r = 1:	Inverse Gamma with c0, C0
  ## r > 1:	Wishart with c0, C0
  ## any r:	Q in {Inverse Gamma, Inverse Wishart} with prQnu (prior Q nu) and prQS
  ## 		We use the Gamma and Wishart and sample the inverse Variance.
  ## where:
  ##      prQnu:	    degrees of freedom for Wishart and rate for Gamma
  ## 		prQS :  	shape for Q
  ##
  ## Select Q0 the prior mean of Q.
  ## Determine prQS from prQS = Q0 * (prQnu - (r + 1)/2). This matches Q0 to the mean
  ## of the Inverse Gamma or the Inverse Wishart distribution and to the mode of Q^{-1}
  ## i.e. the Gamma and Wishart distribution respectively.
  ## Further, variance shrinkage towards the ratio prQS/dfQpr, where dfQpr bounds the
  ## ratio of the variances.

  dfQpr <- 2.5 ## this bounds the ratio of variances to 10 for r = 1
  prQnu <- dfQpr + (r - 1) / 2

  if (K == 1) {
    phi <- 1 ## c0 heterogeneity
  } else {
    ## Tuning of the prior for sigma is done by explained heterogeneity
    ## See p. 192, chapter 6.3.2 Fruewirth-Schnatter (2006)
    ## Rhet
    ## -> 1: 	means very different in relation to variances
    ## -> 0: 	means rather similar in relation to variances
    ## 0 < Rhet < 1 (do not choose 0 nor 1)
    ## SMALL VALUE: 	leads to very informative prior for mu_k
    ## 			close to b0. Should be chosen only in
    ## 			combination with a hierarchical prior
    ## 			on b0.
    ## LARGE VALUE:		leads to a very informative prior for
    ## 			sigma_k close to prQS/prQnu. Should only
    ## 			be chosen in combination with hierarchical
    ## 			prior on prQS.
    Rhet <- 0.5 ## Rhet = 2/3
    phi <- (1 - Rhet)
  }
  prQS <- cov * phi * (prQnu - (r + 1) / 2)
  if (r > 1) {
    detprQS <- log(det(prQS))
  }
  if (hier) {
    if (rich.green) {
      if (r == 1) {
        g0 <- 0.2 ## Richardson and Green. Sampling from Gamma allows
        ## arbitrary g0:
        ## WARNING: seems to cause problems in bayesf
        prQnu <- 2 ## Note that prQnu standard is changed here
      } else if (r == 2) {
        g0 <- 0.3 ## Stephens
        ## WARNING:  seems to cause problems in bayesf
        prQnu <- 3 ## prQnu is changed also in relation from standard
      } else { ## r > 2
        g0 <- 0.5 + (r - 1) / 2
      }
      g0 <- 0.5 + (r - 1) / 2
      G0 <- 100 * g0 / prQnu * solve(cov) ## Stephens
      prQS <- prQnu * cov / 100 ## define starting values for prQS
    } else { ## Bensmail et al.
      g0 <- 0.5 + (r - 1) / 2 ## in general g0 must be a multiple of 0.5 for the
      ## Inverse Wishart (IW) to lead to a proper prior
      G0 <- g0 * solve(prQS) ## match hierarchical and non-hierarchical priors
    }
    if (r > 1) {
      par.sigma <- list(
        c = array(prQnu, dim = c(1, K)),
        C = array(prQS, dim = c(r, r, K)),
        logdetC = array(detprQS, dim = c(1, K)),
        g = g0, G = G0
      )
    } else { ## r == 1
      par.sigma <- list(
        c = array(prQnu, dim = c(1, K)),
        C = array(prQS, dim = c(1, K)),
        g = g0, G = G0
      )
    }
  } else { ## non-hierarchical prior
    if (r > 1) {
      par.sigma <- list(
        c = array(prQnu, dim = c(1, K)),
        C = array(prQS, dim = c(r, r, K)),
        logdetC = array(detprQS, dim = c(1, K))
      )
    } else { ## r == 1
      ## later distinguish between 'sigmauniform' and 'others' ##
      par.sigma <- list(
        c = array(prQnu, dim = c(1, K)),
        C = array(prQS, dim = c(1, K))
      )
    }
  }
  obj@par <- list(mu = par.mu, sigma = par.sigma)
  return(obj)
}

#' Generate default prior for the degrees of freedom
#' 
#' @description 
#' For internal usage only. This function constructs the default priors for 
#' the degrees of freedom.
#' 
#' @param object A `prior` object to be enriched by prior parameters. 
#' @return A `prior` object with specified prior parameters for the prior of 
#'   a Student-t mixture model.
#' @noRd
#' @seealso 
#' * [priordefine()] for the calling function.
".generateDfPrior" <- function(object) {
  ## default prior: independent hierarchical prior following Fernandéz and Steel (1999)
  df.type <- "inhier"
  df.trans <- 1
  df.a0 <- 2
  df.b0 <- 2
  df.mean <- 10
  df.d <- (df.mean - df.trans) * (df.b0 - 1)
  df <- list(
    type = df.type, trans = df.trans,
    a0 = df.a0, b0 = df.b0, d = df.d
  )
  object@par$df <- df
  return(object)
}

#' Generate default prior for the weights of any finite mixture
#' 
#' @description 
#' For internal usage only. This function constructs the default priors for 
#' the weights of any finite mixture.
#' 
#' @details 
#' e_1,...e_K,  1 x K
#' A default with e_i = 4 for all i = 1, ..., K is chosen.
#' 
#' @param object A `prior` object to be enriched by prior parameters. 
#' @param model A `model` object specifying the mixture model.
#' @return A `prior` object with specified prior parameters for the weights of 
#'   any finite mixture model.
#' @noRd
#' 
#' @seealso 
#' * [priordefine()] for the calling function.
".generatePriorWeight" <- function(object, model) {
  K <- model@K
  if (K > 1 && !model@indicfix) {
    e0 <- 4
    object@weight <- matrix(e0, nrow = 1, ncol = K)
  } else { ## K = 1
    object@weight <- matrix()
  }
  return(object)
}

#' Check validity of type of a prior
#' 
#' @description 
#' For internal usage only. This function checks the `type` of a prior. Only 
#' two values are possible: either `"condconjugate"` indicating a conditionally 
#' conjugate prior or `"independent"` for an independence prior.
#' 
#' @param obj A `prior` object with defined `@@type` slot.
#' @return None. Throws an error if the type is wrong.
#' @noRd
#' 
#' @seealso 
#' * [prior()] for the class constructor calling this checking function
#' * [priordefine()] for the advanced class constructor calling this function
".valid.type.Prior" <- function(obj) {
  type.choices <- c("condconjugate", "independent")
  if (!(obj@type %in% type.choices)) {
    stop(paste("Unknown prior 'type'. 'type' must be",
      "'independent' or 'condconjugate'.",
      sep = ""
    ))
  }
  #    if (model.obj@dist == "poisson" && obj@type == "independent") {
  #        warning(paste("Wrong specification of slot 'type' in 'prior' ",
  #                      "object with slot 'dist' in 'model' object set to ",
  #                      "'poisson'. For Poisson mixtures only the prior ",
  #                      "type 'condconjugate' is available.", sep = ""))
  #    }
}

#' Check validity of type of a prior
#' 
#' @description 
#' For internal usage only. This function checks the argument `coef.mat` in 
#' the constructors. This argument is not yet implemented for usage by the user.
#' 
#' @details 
#' The coefficient matrix `coef.mat` for `cond.poisson`
#' distributions with conditional prior must be a lower
#' triangular matrix with ones on its diagonal.
#' Further it must be of type `matrix` or `array` with
#' dimension `KxK`.
#' 
#' @param model.obj A `model` object specifying the mixture model.
#' @param coef.mat A `matrix` containing the coefficients. Must be lower 
#'   triangular.
#' @return None. Throws an error if the coefficient matrix is wrongly specified.
#' @noRd 
#' 
#' @seealso 
#' * [prior()] for the class constructor calling this checking function
#' * [priordefine()] for the advanced class constructor calling this function
".valid.coefmat.Prior" <- function(model.obj, coef.mat) {
  K <- model.obj@K
  if (is.null(coef.mat)) {
    stop("For a conditional Poisson mixture a coefficient matrix
             'coef.mat' has to be provided.")
  } else if (!is.null(coef.mat)) {
    if (!is.matrix(coef.mat) && !is.array(coef.mat)) {
      stop("Argument 'coef.mat' must be of type 'matrix' or 'array'.")
    } else if (nrow(coef.mat) != ncol(coef.mat)) {
      stop("Argument 'coef.mat' must be a quadratic 'matrix' or 'array'.")
    } else if (nrow(coef.mat) != K || ncol(coef.mat) != K) {
      stop("Dimension of argument 'coef.mat' must correspond to number
                 of components 'K' in 'model'.\n")
    } else if (!(all(diag(coef.mat) == 1))) {
      stop("Coefficients on the diagonal of 'coef.mat' must be equal
                 to one.\n")
    }
  }
}

#' Check consistency of `fdata` and `model` object for a prior
#' 
#' @description 
#' For internal usage only. This function checks the consistency of an `fdata` 
#' object and a corresponding `model` object. Consistency is ensured, if the 
#' distribution in slot `@@dist` of a `model` object conforms to the dimension 
#' in slot `@@r` of the `fdata` object.
#' 
#' @param fdata.obj An `fdata` object containing the observations.
#' @param model.obj A `model` object specifying the finite mixture model.
#' @return None. Throws an error if the type is wrong.
#' @noRd
#' 
#' @seealso 
#' * [prior()] for the class constructor calling this checking function
#' * [priordefine()] for the advanced class constructor calling this function 
".valid.fdata.model.Prior" <- function(fdata.obj, model.obj) {
  if (model.obj@dist %in% .get.univ.Model() && fdata.obj@r > 1) {
    stop(paste("Wrong specification of slot 'r' in 'fdata' object. ",
      "Univariate distribution in slot 'dist' of 'model' ",
      "object but dimension in slot 'r' of 'fdata' object ",
      "greater 1.",
      sep = ""
    ))
  } else if (model.obj@dist %in% .get.multiv.Model() && fdata.obj@r < 2) {
    stop(paste("Wrong specification of slot 'r' in 'fdata' object. ",
      "Multivariate distribution in slot 'dist' of 'model' ",
      "object but dimension in slot 'r' of 'fdata' object ",
      "less than two.",
      sep = ""
    ))
  }
}
