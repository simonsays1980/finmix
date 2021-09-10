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

### ================================================================
### The prior class
### ----------------------------------------------------------------
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

### ----------------------------------------------------------------
### prior
### @description    Default constructor.
### @par    weight  an R 'matrix' object containing the prior weights
### @par    par     an R list object containing the hyper parameters
### @par    type    an R 'character' object defining the type of the
###                 prior; possible type are either "independent" or
###                 "condconjugate"
### @par    hier    an R 'logical' object indicating if a hierarchical
###                 prior should be used
### @returns        an S4 object of class 'prior'
### @see    ?prior
### @author Lars SImon Zehnder
### -----------------------------------------------------------------
"prior" <- function(weight = matrix(), par = list(),
                    type = c("independent", "condconjugate"),
                    hier = TRUE) {
  type <- match.arg(type)
  .prior(
    weight = weight, par = par,
    type = type, hier = hier
  )
}
### -----------------------------------------------------------------
### priordefine
### @description    Advanced constructor. Constructs an object from
###                 input parameters. Constructed prior has data-
###                 dependent hyper parameters.
### @par    fdata       an S4 object of class 'fdata'
### @par    model       an S4 object of class 'model'
### @par    coef.mat    not implemented yet
### @par    varargin    an S4 object of class 'prior'
### @return         an S4 object of class 'prior'
### @see    ?fdata, ?model, ?priordefine
### @author Lars Simon Zehnder
### -----------------------------------------------------------------
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

### -----------------------------------------------------------------------
### generaterPrior
### @description    Generates an object of class 'prior' from input
###                 parameters, i.e. it fills all slots with appropriate
###                 values. The object itself is constructed before this
###                 method is called.
### @par    obj         an S4 object of class 'prior'
### @par    fdata       an S4 object of class 'fdata'
### @par    model       an S4 object of class 'model'
### @par    varargin    an S4 object of class 'prior' or 'missing'
### @par    coef.mat    not yet implemented
### @return         a fully specified S4 object of class 'prior'
### @see    .generatePrior<model@dist>
### @author Lars Simon Zehnder
### -----------------------------------------------------------------------
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
setMethod(
  "getWeight", "prior",
  function(object) {
    return(object@weight)
  }
)

setMethod(
  "getPar", "prior",
  function(object) {
    return(object@par)
  }
)

setMethod(
  "getType", "prior",
  function(object) {
    return(object@type)
  }
)

setMethod(
  "getHier", "prior",
  function(object) {
    return(object@hier)
  }
)

## Setters ##
setReplaceMethod(
  "setWeight", "prior",
  function(object, value) {
    object@weight <- value
    validObject(object)
    return(object)
  }
)

setReplaceMethod(
  "setPar", "prior",
  function(object, value) {
    object@par <- value
    validObject(object)
    return(object)
  }
)

setReplaceMethod(
  "setType", "prior",
  function(object, value) {
    object@type <- value
    validObject(object)
    return(object)
  }
)

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

### ==================================================================
### Checking
### ------------------------------------------------------------------

### ------------------------------------------------------------------
### .check.fdata.model.Prior
### @description    Checks objects of classes 'fdata' and 'model' for
###                 validity and consistency.
### @par    fdata.obj   an S4 object of class 'fdata'
### @par    model.obj   an S4 object of class 'model'
### @return         throws an error if any object is not valid or if
###                 the two objects are not consistent among each other
### @see    fdata:::.valid.Fdata, fdata:::.hasY, model:::.valid.Model,
###         .valid.fdata.model.Prior
### @author Lars Simon Zehnder
### -------------------------------------------------------------------
".check.fdata.model.Prior" <- function(fdata.obj, model.obj) {
  .valid.Fdata(fdata.obj)
  hasY(fdata.obj, verbose = TRUE)
  .init.valid.Model(model.obj)
  .valid.fdata.model.Prior(fdata.obj, model.obj)
}

### ------------------------------------------------------------------
### .check.varargin.Prior
### @description    Checks if the variable argument 'varargin' is also
###                 of class 'prior' and is valid. Throws an error if
###                 any condition is not fulfilled.
### @par    obj     any R object passed to the function 'priordefine()'
###                 by the user
### @return         throws an error if 'obj' is not of class 'prior'
### @see    validity
### @author Lars Simon Zehnder
### -------------------------------------------------------------------
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
                      "'gG'.",
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
              "need Beta shape parameters named ",
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

### -----------------------------------------------------------------
### .haspar.exponential.Prior
### -----------------------------------------------------------------
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
        }
      }
    }
  }
}

### -----------------------------------------------------------------
###  obj, model.obj, verbose .generatePriorPoisson
### @description    Generates the hyper parameters for a Poisson
###                 distribution.
### @par    obj     an S4 object of class 'prior'
### @par    fdata.obj   an S4 object of class 'fdata'
### @par    model.obj   an S4 object of class 'model'
### @par    varargin    am S4 object of class 'prior'
### @return         a fully specified 'prior' object for Poisson
###                 models specified by 'model.obj' with data in
###                 'fdata.obj' and predefined slots in 'varargin'
### @details        the type of a data-dependent Poisson prior is
###                 is always conditionally conjugate Gamma with
###                 parameters:     a: shape,   1 x model.obj@K
###                                 b: rate,    1 x model.obj@K
###                 If not otherwise specified in 'varargin' an
###                 hierarchical Gamma distribution is chosen with
###                 parameters:     g: shape,   1 x 1
###                                 G: rate,    1 x 1.
### @see ?prior, ?fdata, ?model
### @author Lars Simon Zehnder
### ------------------------------------------------------------------
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
### ----------------------------------------------------------------
### .generatePriorBinomial
### @description    Generates the hyper parameters for a Binomial
###                 distribution.
### @par    obj         an S4 object of class 'prior'
### @par    model.obj   an S4 object of class 'model'
### @return         a fully specified 'prior' object for Binomial
###                 models specified by 'model.obj'.
### @details        the type of generated Binomial prior is always
###                 conditionally conjugate Beta with parameters:
###                 a:  shape,  1 x model.obj@K
###                 b:  rate,   1 x model.obj@K;
###                 starting values are a = (1, 1), b = (1, 1).
### @see ?prior, ?model
### author Lars Simon Zehnder
### ----------------------------------------------------------------
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

### ----------------------------------------------------------------
### .generatePriorExponential
### @description    Generates the hyper parameters for an Exponential
###                 distribution.
### @param  obj             an S4 object of class 'model'
### @param  fdata.obj       an S4 object of class 'fdata'
### @param  model.obj       an S4 object of class 'model'
### @param  prior.wagner    an R object of class 'logical'
### @return a fully specified 'prior' object for Exponential models
###         specified by 'model.obj' and data specified in'fdata.obj'
### @detail If the identifier 'prior.wagner == TRUE' the prior from
###         Wagner (2007) is taken. In the remaining case a
###         prior is constructed from the analysis of overdispersion
###         in the observations. This prior can also be hierarchical
###         if specified.
### @see    ?priordefine
### @author Lars Simon Zehnder
### ----------------------------------------------------------------
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

### Prior weight: The prior distribution of the weights.
### @Distribution:      Dirichlet
### @Parameters:
###                     e_1,...e_K,  1 x K
### A default with e_i = 4 for all i = 1, ..., K is chosen.
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

### Validity
### Valid type: The prior @type must be one of the two choices
### 'independent' or 'condconjugate' (conditional conjugate).
### For some distribution models only one type of prior exists:
### @Poisson:   'condconjugate'
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

### The coefficient matrix 'coef.mat' for 'cond.poisson'
### distributions with conditional prior must be a lower
### triangular matrix with ones on its diagonal.
### Further it must be of type 'matrix' or 'array' with
### dimension K x K.
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

### -------------------------------------------------------------------------------
### .valid.fdata.model.Prior
### @description    Checks for consistency between the specified model in slot
###                 @dist of the 'model' object and the dimension of variables
###                 @r in the 'fdata' object. Throws and error if no consistency
###                 exists.
### @par    fdata.obj   an S4 object of class 'fdata'
### @par    model.obj   an S4 object of class 'model'
### @return         Throws an error if no consistency is found.
### @see    ?fdata, ?model
### @author Lars Simon Zehnder
### --------------------------------------------------------------------------------
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
