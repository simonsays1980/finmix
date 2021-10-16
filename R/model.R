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

#' Finmix `model` class
#' 
#' @description 
#' This class specifies a finite mixture model. Entities are created from it by 
#' calling its constructor [model()].
#' 
#' @details 
#' A finite mixture model in the ` finmix` package is defined by its number of 
#' components `K`, the component distributions `dist`, the data dimension `r` 
#' and an indicator defining, if the model has fixed or unknown indicators. 
#' Finite mixture models for the following distributions can be constructed: 
#' 
#'  * Poisson,
#'  * Conditional Poisson,
#'  * Exponential,
#'  * Binomial,
#'  * Normal,
#'  * Multivariate Normal,
#'  * Student-t,
#'  * Multivariate Student-t.
#'  
#' Using the constructor [model()] a finite mixture model can be created, the 
#' default being a mixture model of Poisson distributions. 
#' 
#' ## Fully defined finite mixture models
#' A fully defined finite mixture model contains next to the distribution and 
#' the components also weights and parameters. The weights are defined in slot 
#' `weight` and must be of class ` matrix` with as many weights as there are 
#' components in the mixture model (dimension `Kx1`). Parameters are defined in 
#' a ` list` named `par`. The elements of this list depend on the chosen 
#' distribution in slot `dist`: 
#' 
#'  * Poisson: A `matrix` named `lambda` of dimension `Kx1` holding the rate 
#'    parameters.
#'  * Exponential: A `matrix` named `lambda` of dimension `Kx1` holding the rate 
#'    parameters.
#'  * Binomial: A `matrix` of dimension `Kx1` named `p` storing the 
#'    probabilities.
#'  
#' 
#'  
#' 
#' @slot dist A character, defining the distribution family. Possible choices
#' are binomial, exponential, normal, normult, poisson, student, and studmult.
#' @slot r An integer. Defines the vector dimension of a model. Is one for all
#' univariate distributions and larger than one for normult and studmult.
#' @slot K An integer, defining the number of components in the finite mixture.
#' @slot weight A matrix, containing the weights of the finite mixture model. 
#' The matrix must have dimension \code{1 x K} and weights must add to one
#' must all be larger or equal to zero. 
#' @slot par A list containing the parameter vectors for the finite mixture 
#' distribution. The list can contain more than one named parameter vector. 
#' @slot indicmod A character defining the indicator model. So far only 
#' multinomial indicator models are possible. 
#' @slot indicfix A logical. If \code{TRUE} the indicators are given and
#' therefore fixed. 
#' @slot T A matrix containing the repetitions in case of a \code{"binomial"} or 
#'  \code{"poisson"} model.
#' @exportClass model
#' @rdname model-class
#' 
#' @seealso 
#' * [mixturemcmc()] for performing MCMC sampling with a mixture model
#' * [modelmoments()] for compute theoretical moments of a finite mixture model
.model <- setClass("model",
  representation(
    dist = "character",
    r = "integer",
    K = "integer",
    weight = "matrix",
    par = "list",
    indicmod = "character",
    indicfix = "logical",
    T = "matrix"
  ),
  validity = function(object) {
    .init.valid.Model(object)
    ## else: OK ##
    TRUE
  },
  prototype(
    dist = character(),
    r = integer(),
    K = integer(),
    weight = matrix(),
    par = list(),
    indicmod = character(),
    indicfix = logical(),
    T = matrix()
  )
)

#' Constructor for the S4 model class
#' 
#' \code{model} creates a finite mixture model with given parameters. 
#' 
#' This is a constructor that creates a class object and guides the user in 
#' regard to the different parameters needed to define a finite mixture model.
#' 
#' @param dist A character, defining the distribution family. Possible choices
#' are \code{"binomial"}, \code{"exponential"}, \code{"normal"}, 
#' \code{"normult"}, \code{"poisson"}, \code{"student"}, and \code{"studmult"}.
#' @param r An integer. Defines the vector dimension of a model. Is one for all
#' univariate distributions and larger than one for \code{"normult"} and 
#' \code{"studmult"}.
#' @param K An integer, defining the number of components in the finite mixture. 
#' Must be larger or equal to one. 
#' @param weight A matrix, containing the weights of the finite mixture model. 
#' The matrix must have dimension \code{1 x K} and weights must add to one
#' and must all be larger or equal to zero. 
#' @param par A list containing the parameter vectors for the finite mixture 
#' distribution. The list can contain more than one named parameter vector. 
#' Depending on the distribution parameters must be defined in the list as 
#' follows: a \code{K}-dimensional vector of probabilities named \code{"p"} for 
#' a \code{"binomial"} model, a \code{K}-dimensional vector of positive rates 
#' named \code{"lambda"} for an \code{"exponential"} model, 
#' \code{K}-dimensional vectors of means named \code{"mu"} and variances named 
#' \code{sigma} for a \code{"normal"} model, a \code{r x K}-dimensional 
#' matrix of means named \code{"mu"} and a \code{K x r x r} dimensional
#' array of variance-covariance matrices named \code{"sigma"} for a 
#' \code{"normult"} model, a \code{K}-dimensional vector of rates named 
#' \code{"rates"} for a \code{"poisson"} model, \code{K}-dimensional vectors of 
#' means named \code{"mu"}, variances named \code{sigma}, and degrees of freedom
#'  named \code{"df"} for a \code{"student"} model, a 
#' \code{r x K}-dimensional matrix of means named \code{"mu"}, a 
#' \code{K x r x r} dimensional array of variance-covariance matrices 
#' named \code{"sigma"}, and a \code{K}-dimensional vector of degrees of freedom
#'  for a \code{"studmult"} model.
#' @param indicmod A character defining the indicator model used. For now only
#' \code{"multinomial"} is implemented.
#' @param indicfix A logical. If \code{TRUE} the indicators are given and
#' therefore fixed. 
#' @param T A matrix containing the repetitions in case of a \code{"binomial"} or 
#'  \code{"poisson"} model. Must be positive integers.
#' @return An S4 `model` object.
#' @export
#' 
#' @examples 
#' f_model <- model(dist = "poisson", K = 2, par = list(lambda = c(0.17, 0.2)))
#' 
#' @seealso 
#' * [model][model_class] for the class definition
"model" <- function(dist = "poisson", r, K,
                    weight = matrix(), par = list(),
                    indicmod = "multinomial",
                    indicfix = FALSE, T = matrix()) {
  if (missing(K)) {
    K <- .check.K.Model(weight)
  } else {
    K <- as.integer(K)
    if (K == 1 && dist == "cond.poisson") {
      dist <- "poisson"
    }
  }
  if (missing(r)) {
    r <- .check.r.Model(dist)
  } else {
    r <- as.integer(r)
  }
  if (missing(weight) && K > 1) {
    weight <- .check.weight.Model(K)
  } else {
    weight <- as.matrix(weight)
  }
  if (!missing(T)) {
    T <- .check.T.Model(T)
  } else {
    if (dist == "binomial") {
      T <- matrix(as.integer(1))
    }
  }

  .model(
    dist = dist, r = r, K = K, weight = weight,
    par = par, indicmod = indicmod,
    indicfix = indicfix, T = T
  )
}

#' Getter for weights
#'
#' \code{hasWeight} returns the weight matrix. 
#' 
#' @param model An S4 model object. 
#' @param verbose A logical indicating, if the function should give a print out.
#' @return Matrix of weights.
#' @exportMethod hasWeight
#'
#' @examples
#' \dontrun{
#' weight <- hasWeight(model)
#' }
#' @rdname model_class
setMethod(
  "hasWeight", "model",
  function(object, verbose = FALSE) {
    if (!all(is.na(object@weight))) {
      if (ncol(object@weight) == object@K) {
        return(TRUE)
      } else {
        if (verbose) {
          stop(paste("Wrong dimension of ",
            "slot 'weight' of ",
            "'model' object.",
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
        stop(paste("Slot 'weight' of 'model' ",
          "object is empty.",
          sep = ""
        ))
      } else {
        return(FALSE)
      }
    }
  }
)

#' Checks for repetitions.
#' 
#' \code{hasT} chwecks if the model object possesses repetitions.
#' 
#' @param model An S4 model object.
#' @param verbose A logical indicating if the function should give a print out. 
#' @return A logical. \code{TRUE} if repetitions are existent in the model. If 
#' values of slot \code{T} are \code{NA} it returns \code{FALSE}.
#' @exportMethod hasT
#' 
#' @examples
#' \dontrun{
#' if(hasT(model)) {cat('Has repetitions.')}
#' }
#' 
#' @seealso \code{model}
setMethod(
  "hasT", "model",
  function(object, verbose = FALSE) {
    if (!all(is.na(object@T))) {
      return(TRUE)
    } else {
      if (verbose) {
        stop(paste("Slot 'T' of 'model' ",
          "object is empty.",
          sep = ""
        ))
      } else {
        return(FALSE)
      }
    }
  }
)

#' Checks for parameters.
#' 
#' \code{hasPar} checks if the model has parameters defined. 
#' 
#' @param model An S4 model object.
#' @param verbose A logical indicating, if the function should give a print out. 
#' @return A matrix with repetitions. Can be empty, if no repetitions are set.
#' @exportMethod hasPar
#' 
#' @examples 
#' \dontrun{
#' if(hasPar(model)) {simulate(model)}
#' }
#' 
#' @seealso \code{model}
setMethod(
  "hasPar", "model",
  function(object, verbose = FALSE) {
    .haspar.Model(object, verbose)
  }
)

#' Simulates data from a model. 
#' 
#' `simulate()` simulates values for a specified mixture model in an 
#' S4 `model` object.
#' 
#' @param model An S4 model object with specified parameters and weights.
#' @param N An integer specifying the number of values to be simulated. 
#' @param varargin An S4 fdata object with specified variable dimensions, `r` 
#'   and repetitions `T`. 
#' @param seed An integer specifying the seed for the RNG. 
#' @return An S4 fdata object holding the simulated values.
#' @exportMethod simulate
#' @keywords internal
#' 
#' @examples 
#' \dontrun{
#' f_data <- simulate(model, 100)
#' }
#' 
#' @seealso 
#' * [model-class] for the class definition
#' * [fdata-class] for the class defining `finmix` data objects
setMethod(
  "simulate", "model",
  function(model, N = 100, varargin, seed = 0) {
    ## TODO: Check model for parameters. Check varargin for dimension. Check
    ##      model and varargin for consistency.
    if (!missing(seed)) {
      set.seed(seed)
    } ## Implemented maybe finmixOptions with a state variable seed
    if (!hasWeight(model)) {
      model@weight <- matrix(1 / model@K, nrow = 1, ncol = model@K)
    }
    ## Start simulating the allocations
    S <- .simulate.indicators.Model(model, N)
    if (missing(varargin)) {
      varargin <- fdata(
        r = model@r, T = matrix(1, nrow = N),
        exp = matrix(1, nrow = N), S = S
      )
    } else {
      varargin@S <- S
    }
    if (hasPar(model, verbose = TRUE)) {
      .simulate.data.Model(model, N, varargin)
    }
  }
)

#' Plots a model.
#' 
#' \code{plot} plots the density or probabilities of a fully specified mixture 
#' model.
#' 
#' @param x An S4 model object. Must have specified parameters and weights.
#' @param y Unused.
#' @param dev A logical indicating, if the plot should be shown in a graphical 
#' device. Set to \code{FALSE}, if plotted to a file. 
#' @return Density or barplot of the S4 model object. 
#' @exportMethod plot
#' 
#' @examples \dontrun{
#' plot(f_model)
#' }
#' 
#' @seealso \code{model}
setMethod(
  "plot", "model",
  function(x, y, dev = TRUE, ...) {
    dist <- x@dist
    if (dist == "normal") {
      .plot.Normal.Model(x, dev, ...)
    } else if (dist == "normult") {
      .plot.Normult.Model(x, dev, ...)
    } else if (dist == "exponential") {
      .plot.Exponential.Model(x, dev, ...)
    } else if (dist == "student") {
      .plot.Student.Model(x, dev, ...)
    } else if (dist == "studmult") {
      .plot.Studmult.Model(x, dev, ...)
    } else if (dist %in% c("poisson", "cond.poisson")) {
      .plot.Poisson.Model(x, dev, ...)
    } else if (dist == "binomial") {
      if (abs(max(x@T) - min(x@T)) > 1e-6) {
        stop("Plotting a binomial distribution with varying
                           repetitions in slot 'T' is not possible.")
      }
      .plot.Binomial.Model(x, dev, ...)
    }
  }
)

#' Plots point process.
#' 
#' \code{plotPointProc} plots the point process of an S4 model object that 
#' defines a finite mixture model. Only available for Poisson mixtures so far.
#' 
#' @param x An S4 model object with defined parameters and weigths. 
#' @param y Unused.
#' @param dev A logical indicating, if the plot should be shown in a graphical 
#' device. Set to \code{FALSE}, if plotted to a file. 
#' @return A scatter plot of weighted parameters.  
#' @exportMethod plotPointProc
#' 
#' @examples 
#' \dontrun{
#' plotPointProc(f_model)
#' }
#' 
#' @seealso \code{model}
setMethod(
  "plotPointProc", signature(
    x = "model",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    hasPar(x, verbose = TRUE)
    hasWeight(x, verbose = TRUE)
    if (x@dist == "poisson") {
      .plotpointproc.Poisson(x, dev)
    }
  }
)

## Marginal Mixture ##
#' Returns the marginal distribution. 
#' 
#' \code{mixturemar} returns the marginal distribution of a multivariate 
#' mixture distribution. This can only be applied on S4 model objects with 
#' \code{dist="normult"} or \code{dist="studmult"}. 
#' 
#' @param object An S4 model object with a multivariate distribution.
#' @param J An integer specifying the dimension for which the marginal 
#' distribution should be returned.
#' @return An S4 model object with the marginal distribution for dimension 
#' \code{J}.
#' @exportMethod mixturemar
#' 
#' @examples
#' \dontrun{
#' mar_model <- mixturemar(f_model, 1)
#' }
#' 
#' @seealso \code{model}
setMethod(
  "mixturemar", "model",
  function(object, J) {
    .mixturemar.Model(object, J)
  }
)

#' Shows the model.
#' 
#' \code{show} prints model information to the console. 
#' 
#' @param object An S4 model object. 
#' @return A print out of model information about all slots. 
#' @exportMethod show
#' 
#' @examples 
#' \dontrun{
#' show(f_model)
#' }
#' 
#' @seealso \code{model}
setMethod(
  "show", "model",
  function(object) {
    cat("Object 'model'\n")
    cat("     class       :", class(object), "\n")
    cat("     dist        :", object@dist, "\n")
    cat("     r           :", object@r, "\n")
    cat("     K           :", object@K, "\n")
    if (hasPar(object)) {
      cat(
        "     par         : List of",
        length(object@par), "\n"
      )
    }
    if (!object@indicfix) {
      cat(
        "     weight      :",
        paste(dim(object@weight), collapse = "x"),
        "\n"
      )
    }
    cat("     indicmod    :", object@indicmod, "\n")
    cat("     indicfix    :", object@indicfix, "\n")
    if (object@dist == "binomial" && !all(is.na(object@T))) {
      cat(
        "     T           :",
        paste(dim(object@T), collapse = "x"), "\n"
      )
    }
  }
)

## Getters ##
#' @name model_class
#' @exportMethod getDist
setMethod(
  "getDist", "model",
  function(object) {
    return(object@dist)
  }
)

#' @name model_class
#' @exportMethod getR
setMethod(
  "getR", "model",
  function(object) {
    return(object@r)
  }
)

#' @name model_class
#' @exportMethod getK
setMethod(
  "getK", "model",
  function(object) {
    return(object@K)
  }
)

#' @name model_class
#' @exportMethod getWeight
setMethod(
  "getWeight", "model",
  function(object) {
    return(object@weight)
  }
)

#' @name model_class
#' @exportMethod getPar
setMethod(
  "getPar", "model",
  function(object) {
    return(object@par)
  }
)

#' @name model_class
#' @exportMethod getIndicmod
setMethod(
  "getIndicmod", "model",
  function(object) {
    return(object@indicmod)
  }
)

#' @name model_class
#' @exportMethod getIndicfix
setMethod(
  "getIndicfix", "model",
  function(object) {
    return(object@indicfix)
  }
)

#' @name model_class
#' @exportMethod getT
setMethod(
  "getT", "model",
  function(object) {
    return(object@T)
  }
)

## Setters ##
#' @name model_class
#' @exportMethod setDist<-
setReplaceMethod(
  "setDist", "model",
  function(object, value) {
    object@dist <- value
    .valid.dist.Model(object)
    return(object)
  }
)

#' @name model_class
#' @exportMethod setR<-
setReplaceMethod(
  "setR", "model",
  function(object, value) {
    object@r <- as.integer(value)
    validObject(object)
    return(object)
  }
)

#' @name model_class
#' @exportMethod setK<-
setReplaceMethod(
  "setK", "model",
  function(object, value) {
    object@K <- as.integer(value)
    .valid.K.Model(object)
    if (object@K > 1) {
      object@weight <- .check.weight.Model(object@K)
    } else {
      weight <- matrix()
      storage.mode(weight) <- "numeric"
      object@weight <- weight
    }
    return(object)
  }
)

#' @name model_class
#' @exportMethod setWeight<-
setReplaceMethod(
  "setWeight", "model",
  function(object, value) {
    object@weight <- as.matrix(value)
    object@K <- ncol(object@weight)
    .valid.weight.Model(object)
    return(object)
  }
)

#' @name model_class
#' @exportMethod setPar<-
setReplaceMethod(
  "setPar", "model",
  function(object, value) {
    object@par <- value
    .valid.par.Model(object)
    return(object)
  }
)

#' @name model_class
#' @exportMethod setIndicmod<-
setReplaceMethod(
  "setIndicmod", "model",
  function(object, value) {
    object@indicmod <- value
    return(object)
  }
)

#' @name model_class
#' @exportMethod setIndicfix<-
setReplaceMethod(
  "setIndicfix", "model",
  function(object, value) {
    object@indicfix <- value
    return(object)
  }
)

#' @name model_class
#' @exportMethod setT<-
setReplaceMethod(
  "setT", "model",
  function(object, value) {
    object@T <- matrix(value)
    .valid.T.Model(object)
    return(object)
  }
)

### Private functions
### These functions are not exported

### Checking.
### Checking is used for in the constructor.
### Arguments for the slots are checked for validity and
### if missing are given by default values. Altogether the
### constructor tries to construct a fully specified model
### object with consistent slots.

### Check K: If weights are provided by the user, the number
### of components is set to the number of columns of the weights.
### If argument 'weight' is missing from the call, the number of
### components is assumed to be one.

#' @noRd
".check.K.Model" <- function(weight) {
  if (!all(is.na(weight))) {
    return(NCOL(weight))
  } else {
    return(as.integer(1))
  }
}

### Check r: The dimension of the model is determined in regard to
### the defined distribution in argument 'dist' (if missing the
### default is 'poisson'). For univariate distributions it is set
### to one and for multivariate distribution as a default to two.
#' @noRd
".check.r.Model" <- function(dist) {
  univ <- .get.univ.Model()
  multiv <- .get.multiv.Model()
  if (dist %in% univ) {
    return(as.integer(1))
  } else if (dist %in% multiv) {
    return(as.integer(2))
  } else {
    stop(paste("Unknown distribution in slot ",
      "'dist' of 'model' object.",
      sep = ""
    ))
  }
}

### Check weight: If argument 'weight' is missing from the call
### equally balanced weights are given as a default.
#' @noRd
".check.weight.Model" <- function(K) {
  weight <- matrix(1 / K, nrow = 1, ncol = K)
  return(weight)
}

### Check T: If repetitions are given they are checked in regard
### to validity. In case of non-numeric objects an error is thrown.
### In case of objects of type 'numeric' it is implicitly converted
### to type 'integer'.
#' @noRd
".check.T.Model" <- function(T) {
  if (!all(is.na(T))) {
    if (!is.numeric(T)) {
      stop(paste("Wrong specification of slot 'T' in ",
        "'model' object. Repetitions must be of ",
        "type 'integer'.",
        sep = ""
      ))
    } else {
      storage.mode(T) <- "integer"
      return(T)
    }
  }
}

### Marginal model
#' @noRd
".mixturemar.Model" <- function(obj, J) {
  if (obj@dist == "normult") {
    .mixturemar.normult.Model(obj, J)
  } else if (obj@dist == "studmult") {
    .mixturemar.studmult.Model(obj, J)
  } else {
    stop("A marginal distribution can only be obtained from 
         multivariate distributions.")
  }
}

#' @noRd
".mixturemar.normult.Model" <- function(obj, J) {
  dist <- ifelse(length(J) == 1, "normal", "normult")
  r <- length(J)
  K <- obj@K
  weight <- obj@weight
  mu <- obj@par$mu[J, ]
  sigma <- obj@par$sigma[J, J, ]
  par <- list(mu = mu, sigma = sigma)
  indicmod <- "multinomial"
  indicfix <- TRUE
  margin.model <- .model(
    dist = dist, r = r, K = K,
    weight = weight, par = par,
    indicmod = indicmod,
    indicfix = indicfix
  )
  validObject(margin.model)
  return(margin.model)
}

#' @noRd
".mixturemar.studmult.Model" <- function(obj, J) {
  dist <- ifelse(length(J) == 1, "student", "studmult")
  r <- length(J)
  K <- obj@K
  weight <- obj@weight
  mu <- obj@par$mu[J, ]
  sigma <- obj@par$sigma[J, J, ]
  df <- obj@par$df
  par <- list(mu = mu, sigma = sigma, df = df)
  indicmod <- "multinomial"
  indicfix <- TRUE
  margin.model <- .model(
    dist = dist, r = r, K = K,
    weight = weight, par = par,
    indicmod = indicmod,
    indicfix = indicfix
  )
  validObject(margin.model)
  return(margin.model)
}

### ==============================================================
### Simulate
### --------------------------------------------------------------

### --------------------------------------------------------------
### .simulate.indicators.Model
### @description    Simulates the indicators.
### @par    obj an S4 object of class 'model'
### @par    N   an R 'integer' object
### @return         an R 'matrix' object with N simulated indi-
###                 cators.
### @details        indicators are simulated via the slot @weight
###                 the 'model' object
### @see    ?simulate
### @author Lars Simon Zehnder
### --------------------------------------------------------------

### TODO: Implement C++ function.
#' @noRd
".simulate.indicators.Model" <- function(obj, N) {
  K <- obj@K
  if (K == 1) {
    S <- matrix(as.integer(1), nrow = N, ncol = K)
  } else {
    ## if (model@indicmod = "") -> "Multinomial"
    ## if Markov else
    if (obj@indicmod == "multinomial") {
      rnd <- runif(N)
      rnd <- matrix(rnd, nrow = N, ncol = K)
      weightm <- matrix(obj@weight,
        nrow = N, ncol = K,
        byrow = TRUE
      )
      S <- apply((t(apply(weightm, 1, cumsum)) < rnd), 1, sum) + 1
      S <- matrix(S, nrow = N)
      storage.mode(S) <- "integer"
    }
  }
  return(S)
}

### --------------------------------------------------------------------
### .simulate.data.Model
### @description    Simulates the simulation functions for a specific model.
### @par    obj         an S4 'model' object
### @par    N           an R 'integer' object; number of simulated values
### @par    fdata.obj   an S4 'fdata' object
### @return         an S4 object of class 'fdata' with simulated values
### @see    ?fdata, ?simulate
### @author Lars Simon Zehnder
### ---------------------------------------------------------------------
#' @noRd
".simulate.data.Model" <- function(obj, N, fdata.obj) {
  dist <- obj@dist
  if (dist == "poisson" || dist == "cond.poisson") {
    .simulate.data.poisson.Model(obj, N, fdata.obj)
  } else if (dist == "binomial") {
    .simulate.data.binomial.Model(obj, N, fdata.obj)
  } else if (dist == "exponential") {
    .simulate.data.exponential.Model(obj, N, fdata.obj)
  } else if (dist == "normal") {
    .simulate.data.normal.Model(obj, N, fdata.obj)
  } else if (dist == "student") {
    .simulate.data.student.Model(obj, N, fdata.obj)
  } else if (dist == "normult") {
    .simulate.data.normult.Model(obj, N, fdata.obj)
  }
}

### ---------------------------------------------------------------------
### .simulate.data.poisson.Model
### @description    Simulates values from a Poisson mixture using pre-
###                 specified model and indicators
### @par    obj         an S4 object of class 'model'
### @par    N           an R 'integer' object; number of simulated values
### @par    fdata.obj   an S4 object of class 'fdata'
### @return         an S4 object of class 'fdata' with simulated values
### @see    ?simulate, model:::.simulate.data.Model, ?rpois
### @author Lars Simon Zehnder
### ---------------------------------------------------------------------
#' Simulate data from a Poisson finite mixture model
#' 
#' @description 
#' Simulates values from a Poisson mixture using pre-specified model and 
#' indicators.
#' 
#' @param obj A `model` object specifying the finite mixture model. 
#' @param N An integer specifying the sample size.
#' @param fdata.obj An `fdata` object to store the simulated data.
#' @return An `fdata` object with simulated data.
#' @importFrom stats rpois
#' @noRd
".simulate.data.poisson.Model" <- function(obj, N, fdata.obj) {
  fdata.obj@type <- "discrete"
  fdata.obj@sim <- TRUE
  fdata.obj@y <- matrix(rpois(N, fdata.obj@exp * obj@par$lambda[fdata.obj@S]))
  return(fdata.obj)
}

#' Simulate data from Binomial mixture model
#' 
#' @description 
#' Simulates values from a Binomial mixture using pre-specified model and 
#' indicators
#' @param obj An `model` object specifying the mixture model.
#' @param N An integer specifying the size of the simulated sample.
#' @param fdata.obj An `fdata` object to store the simulated sample. If the 
#'   `fdata` object contains repetitions in slot `@@T`, the repetitions are 
#'   used in sampling.
#' @return An `fdata` object containing the simulated values.
#' @importFrom stats rbinom
#' @noRd
#' 
#' @seealso 
#' [simulate()][model_class] for the calling function
".simulate.data.binomial.Model" <- function(obj, N, fdata.obj) {
  if (!hasT(fdata.obj)) {
    fdata.obj@T <- as.matrix(1)
  }
  fdata.obj@type <- "discrete"
  fdata.obj@sim <- TRUE
  fdata.obj@y <- matrix(rbinom(N, fdata.obj@T, obj@par$p[fdata.obj@S]))
  return(fdata.obj)
}

#' Simulate data from exponential mixture model
#' 
#' @description 
#' Simulates values from a exponential mixture using pre-specified model and 
#' indicators
#' @param obj An `model` object specifying the mixture model.
#' @param N An integer specifying the size of the simulated sample.
#' @param fdata.obj An `fdata` object to store the simulated sample. 
#' @return An `fdata` object containing the simulated values.
#' @importFrom stats rexp
#' @noRd
#' 
#' @seealso 
#' [simulate()][model_class] for the calling function
".simulate.data.exponential.Model" <- function(obj, N, fdata.obj) {
  fdata.obj@type <- "continuous"
  fdata.obj@sim <- TRUE
  fdata.obj@y <- matrix(rexp(N, obj@par$lambda[fdata.obj@S]))
  return(fdata.obj)
}

#' Simulate data from Normal mixture model
#' 
#' @description 
#' Simulates values from a Normal mixture using pre-specified model and 
#' indicators
#' @param obj An `model` object specifying the mixture model.
#' @param N An integer specifying the size of the simulated sample.
#' @param fdata.obj An `fdata` object to store the simulated sample. 
#' @return An `fdata` object containing the simulated values.
#' @noRd
#' 
#' @seealso 
#' [simulate()][model_class] for the calling function
".simulate.data.normal.Model" <- function(obj, N, fdata.obj) {
  fdata.obj@type <- "continuous"
  fdata.obj@sim <- TRUE
  fdata.obj@y <- matrix(rnorm(
    N, obj@par$mu[fdata.obj@S],
    obj@par$sigma[fdata.obj@S]
  ))
  return(fdata.obj)
}

#' Simulate data from Student-t mixture model
#' 
#' @description 
#' Simulates values from a Student-t mixture using pre-specified model and 
#' indicators
#' @param obj An `model` object specifying the mixture model.
#' @param N An integer specifying the size of the simulated sample.
#' @param fdata.obj An `fdata` object to store the simulated sample. If the 
#'   `fdata` object contains repetitions in slot `@@T`, the repetitions are 
#'   used in sampling.
#' @return An `fdata` object containing the simulated values.
#' @importFrom stats rgamma
#' @noRd
#' 
#' @seealso 
#' [simulate()][model_class] for the calling function
".simulate.data.student.Model" <- function(obj, N, fdata.obj) {
  fdata.obj@type <- "continuous"
  fdata.obj@sim <- TRUE
  omega <- rgamma(N, obj@par$df[fdata.obj@S] / 2,
    rate = 2 / obj@par$df[fdata.obj@S]
  )
  fdata.obj@y <- as.matrix(obj@par$mu[fdata.obj@S] +
    sqrt(obj@par$sigma[fdata.obj@S] / omega) *
      rnorm(N, 0.0, 1.0))
  return(fdata.obj)
}

#' Simulate data from a multivariate Normal mixture model
#' 
#' @description 
#' Simulates values from a multivariate Normal mixture using pre-specified 
#' model and indicators
#' @param obj An `model` object specifying the mixture model.
#' @param N An integer specifying the size of the simulated sample.
#' @param fdata.obj An `fdata` object to store the simulated sample. If the 
#'   `fdata` object contains repetitions in slot `@@T`, the repetitions are 
#'   used in sampling.
#' @return An `fdata` object containing the simulated values.
#' @importFrom mvtnorm rmvnorm
#' @noRd
#' 
#' @seealso 
#' [simulate()][model_class] for the calling function
".simulate.data.normult.Model" <- function(obj, N, fdata.obj) {
  fdata.obj@type <- "continuous"
  fdata.obj@sim <- TRUE
  fdata.obj@y <- matrix(numeric(), nrow = N, ncol = obj@r)
  fdata.obj@r <- obj@r
  for (i in 1:N) {
    fdata.obj@y[i, ] <- rmvnorm(1,
      mean = obj@par$mu[, fdata.obj@S[i]],
      sigma = obj@par$sigma[, , fdata.obj@S[i]],
      method = "chol"
    )
  }
  return(fdata.obj)
}

### Plotting
### Plot Poisson models: Poisson models are discrete
### models and a barplot is used.
### The range for the x-axis is determined via the
### quantiles of the largest and smallest Poisson model
### in the mixture.
#' @importFrom stats qpois dpois
#' @importFrom grDevices axisTicks
#' @noRd 
".plot.Poisson.Model" <- function(model.obj, dev, ...) {
  if (.check.grDevice() && dev) {
    dev.new(title = "Model plot")
  }
  lambda <- model.obj@par$lambda
  weight <- model.obj@weight
  xlim.up <- qpois(.9999, lambda = max(lambda))
  xlim.low <- qpois(.0001, lambda = min(lambda))
  x.grid <- seq(xlim.low, xlim.up, by = 1)
  y.grid <- sapply(x.grid, dpois, lambda = lambda)
  y.grid <- weight %*% y.grid
  main.title <- paste("Poisson Mixture K = ",
    model.obj@K,
    sep = ""
  )
  label.grid <- axisTicks(c(xlim.low, xlim.up),
    log = FALSE,
    nint = 10
  )
  bp <- barplot(y.grid,
    main = main.title, axes = F,
    col = "gray65", border = "gray65", ...
  )
  axis(side = 2, cex = .7, cex.axis = .7)
  axis(
    side = 1, tick = FALSE, at = bp[which(x.grid %in% label.grid)],
    labels = which(x.grid %in% label.grid), cex.axis = .7
  )
  mtext(side = 1, "x", cex = .7, cex.axis = .7, line = 3)
  mtext(side = 2, "P(x)", cex = .7, cex.axis = .7, line = 3)
}

### Plot Binomial models: Binomial models are discrete
### models and line model is used.
### The grid for the x-axis is determined by taking
### the
#' @importFrom stats dbinom
#' @noRd
".plot.Binomial.Model" <- function(model.obj, dev, ...) {
  if (.check.grDevice() && dev) {
    dev.new(title = "Model plot")
  }
  n <- model.obj@T[1]
  p <- model.obj@par$p
  weight <- model.obj@weight
  xlim <- max(n, na.rm = TRUE)
  x.grid <- seq(0, xlim, by = 1)
  y.grid <- sapply(x.grid, dbinom, size = n, p = p)
  y.grid <- weight %*% y.grid
  main.title <- paste("Binomial Mixture K = ",
    model.obj@K,
    sep = ""
  )
  plot(x.grid, y.grid,
    main = main.title, type = "h",
    xlab = "x", ylab = "P(x)", ...
  )
  points(x.grid, y.grid, pch = 20)
}

#' @importFrom stats qexp dexp
#' @noRd
".plot.Exponential.Model" <- function(model.obj, dev, ...) {
  if (.check.grDevice() && dev) {
    dev.new(title = "Model plot")
  }
  lambda <- model.obj@par$lambda
  weight <- model.obj@weight
  min.lambda <- min(lambda, na.rm = TRUE)
  xlim <- qexp(.9999, rate = min.lambda)
  x.grid <- seq(0, ceiling(xlim),
    length =
      as.integer(100 * lambda^(-2))
  )
  y.grid <- sapply(x.grid, dexp, rate = lambda)
  y.grid <- weight %*% y.grid
  main.title <- paste("Exponential Mixture K = ",
    model.obj@K,
    sep = ""
  )
  plot(x.grid, y.grid,
    main = main.title, type = "l",
    xlab = "x", ylab = "P(x)", ...
  )
}

#' @importFrom stats qt dt
#' @noRd
".plot.Student.Model" <- function(model.obj, dev, ...) {
  if (.check.grDevice() && dev) {
    dev.new(title = "Model plot")
  }
  mu <- model.obj@par$mu
  sigma <- model.obj@par$sigma
  df <- model.obj@par$df
  weight <- model.obj@weight
  max.mu <- max(mu, na.rm = TRUE)
  max.sigma <- max(sigma, na.rm = TRUE)
  min.df <- min(df, na.rm = TRUE)
  xlim <- max.mu + max.sigma * qt(.9999, min.df)
  x.grid <- seq(-xlim, xlim, length = 1000) + max.mu
  y.grid <- sapply(x.grid, "-", mu)
  y.grid <- apply(y.grid, 2, "/", sigma)
  y.grid <- apply(y.grid, 2, dt, df = df)
  y.grid <- apply(y.grid, 2, "/", sqrt(sigma))
  y.grid <- t(weight %*% y.grid)
  main.title <- paste("Student-t Mixture K = ",
    model.obj@K,
    sep = ""
  )
  plot(x.grid, y.grid,
    main = main.title, type = "l",
    xlab = "x", ylab = "P(x)", ...
  )
}

#' @importFrom stats qnorm dnorm
#' @noRd
".plot.Normal.Model" <- function(model.obj, dev, ...) {
  if (.check.grDevice() && dev) {
    dev.new(title = "Model Plot")
  }
  mu <- model.obj@par$mu
  sigma <- model.obj@par$sigma
  weight <- model.obj@weight
  max.mu <- max(mu, na.rm = TRUE)
  max.sigma <- max(mu, na.rm = TRUE)
  xlim <- qnorm(.9999,
    mean = max.mu,
    sd = max.sigma
  )
  x.grid <- seq(-xlim, xlim, length = 1000) + max.mu
  y.grid <- sapply(x.grid, dnorm,
    mean = mu,
    sd = sigma
  )
  y.grid <- weight %*% y.grid
  main.title <- paste("Normal Mixture K = ",
    model.obj@K,
    sep = ""
  )
  plot(x.grid, y.grid,
    main = main.title, type = "l",
    xlab = "x", ylab = "P(x)", ...
  )
}

#' @noRd
".plot.Normult.Model" <- function(model.obj, dev, ...) {
  K <- model.obj@K
  r <- model.obj@r
  if (r == 2) {
    if (.check.grDevice() && dev) {
      dev.new(title = "Model: Perspective plot")
    }
    xyz.grid <- .generate.Grid.Normal(model.obj)
    main.title <- paste("Multivariate Normal Mixture K = ",
      K,
      sep = ""
    )
    persp(xyz.grid$x, xyz.grid$y, xyz.grid$z,
      col = "gray65",
      border = "gray47", theta = 55, phi = 30, expand = 0.5,
      lphi = 180, ltheta = 90, r = 40, d = 0.1,
      ticktype = "detailed", zlab = "P(x)", xlab = "r = 1",
      ylab = "r = 2", cex = 0.7, cex.lab = 0.7, cex.axis = 0.7
    )
  } else if (r > 2 && r < 6) {
    if (.check.grDevice() && dev) {
      dev.new(title = "Model: Contour plots")
    }
    if (r == 3) {
      par(
        mfrow = c(1, r), mar = c(2, 2, 2, 2),
        oma = c(4, 5, 1, 5)
      )
    } else if (r == 4) {
      par(
        mfrow = c(2, 3), mar = c(2, 2, 2, 2),
        oma = c(4, 5, 1, 5)
      )
    } else {
      par(
        mfrow = c(2, 5), mar = c(2, 2, 2, 2),
        oma = c(4, 5, 1, 5)
      )
    }
    for (i in seq(1, r - 1)) {
      for (j in seq(1, r)) {
        marmodel <- mixturemar(model.obj, J = c(i, j))
        xyz.grid <- .generate.Grid.Normal(marmodel)
        contour(xyz.grid$x, xyz.grid$y, xyz.grid$z,
          col = "gray47", cex = 0.7, cex.axis = 0.7,
          xlab = paste("r = ", i, sep = ""),
          ylab = paste("r = ", j, sep = "")
        )
      }
    }
  } else {
    stop("Method 'plot' for 'model' objects is not implemented for
             model dimensions of r > 5.")
  }
}

#' @noRd
".plot.Studmult.Model" <- function(model.obj, dev, ...) {
  K <- model.obj@K
  r <- model.obj@r
  if (r == 2) {
    if (.check.grDevice() && dev) {
      dev.new(title = "Model: Perspective plot")
    }
    xyz.grid <- .generate.Grid.Student(model.obj)
    main.title <- paste("Multivariate Student-t Mixture K = ",
      K,
      sep = ""
    )
    persp(xyz.grid$x, xyz.grid$y, xyz.grid$z,
      col = "gray65",
      border = "gray47", theta = 55, phi = 30, expand = 0.5,
      lphi = 180, ltheta = 90, r = 40, d = 0.1,
      ticktype = "detailed", zlab = "P(x)", xlab = "r = 1",
      ylab = "r = 2", cex = 0.7, cex.lab = 0.7, cex.axis = 0.7
    )
  } else if (r > 2 && r < 6) {
    if (.check.grDevice() && dev) {
      dev.new(title = "Model: Contour plots")
    }
    if (r == 3) {
      par(
        mfrow = c(1, r), mar = c(2, 2, 2, 2),
        oma = c(4, 5, 1, 5)
      )
    } else if (r == 4) {
      par(
        mfrow = c(2, 3), mar = c(2, 2, 2, 2),
        oma = c(4, 5, 1, 5)
      )
    } else {
      par(
        mfrow = c(2, 5), mar = c(2, 2, 2, 2),
        oma = c(4, 5, 1, 5)
      )
    }
    for (i in seq(1, r - 1)) {
      for (j in seq(1, r)) {
        marmodel <- mixturemar(model.obj, J = c(i, j))
        xyz.grid <- .generate.Grid.Student(marmodel)
        contour(xyz.grid$x, xyz.grid$y, xyz.grid$z,
          col = "gray47", cex = 0.7, cex.axis = 0.7,
          xlab = paste("r = ", i, sep = ""),
          ylab = paste("r = ", j, sep = "")
        )
      }
    }
  } else {
    stop("Method 'plot' for 'model' objects is not implemented for
             model dimensions of r > 5.")
  }
}

#' @importFrom mvtnorm qmvnorm dmvnorm
#' @noRd
".generate.Grid.Normal" <- function(model.obj) {
  K <- model.obj@k
  mu <- model.obj@par$mu
  sigma <- model.obj@par$sigma
  weight <- model.obj@weight
  func <- function(s, t) {
    value <- 0
    for (k in seq(1, K)) {
      value <- value + weight[k] *
        dmvnorm(cbind(s, t),
          mean = mu[, k],
          sigma = sigma[, , k]
        )
    }
  }
  mu.norm <- apply(mu, 2, function(x) sqrt(sum(x^2)))
  max.mu.index <- tail(sort(mu.norm, index = TRUE)$ix, 1)
  max.mu <- mu[, max.mu.index]
  sigma.det <- apply(sigma, 3, det)
  max.sigma.index <- tail(sort(sigma.det, index = TRUE)$ix, 1)
  max.sigma <- sigma[, , max.sigma.index]
  xylim <- qmvnorm(.9999,
    mean = max.mu,
    sigma = max.sigma
  )$quantile
  x.grid <- seq(-xylim, xylim, length = 100)
  xy.grid <- cbind(x.grid, x.grid)
  xy.grid <- t(apply(xy.grid, 1, "+", max.mu))
  z.grid <- outer(xy.grid[, 1], xy.grid[, 2], func)
  grid.list <- list(
    x = xy.grid[, 1], y = xy.grid[, 2],
    z = z.grid
  )
  return(grid.list)
}

#' @importFrom mvtnorm qmvt dmvt
#' @noRd
".generate.Grid.Student" <- function(model.obj) {
  K <- model.obj@K
  mu <- model.obj@par$mu
  sigma <- model.obj@par$sigma
  df <- model.obj@par$df
  weight <- model.obj@weight
  func <- function(s, t) {
    value <- 0
    for (k in seq(1, K)) {
      value <- value + weight[k] *
        dmvt(cbind(s, t),
          delta = mu[, k],
          sigma = sigma[, , k], df = df[k]
        )
    }
  }
  mu.norm <- apply(mu, 2, function(x) sqrt(sum(x^2)))
  max.mu.index <- tail(sort(mu.norm, index = TRUE)$ix, 1)
  max.mu <- mu[, max.mu.index]
  sigma.det <- apply(sigma, 3, det)
  max.sigma.index <- tail(sort(sigma.det, index = TRUE)$ix, 1)
  max.sigma <- sigma[, , max.sigma.index]
  min.df <- min(df, na.rm = TRUE)
  xylim <- qmvt(.9999,
    delta = max.mu,
    sigma = max.sigma, df = min.df
  )$quantile
  x.grid <- seq(-xylim, xylim, length = 100)
  xy.grid <- cbind(x.grid, x.grid)
  xy.grid <- t(apply(xy.grid, 1, "+", max.mu))
  z.grid <- outer(xy.grid[, 1], xy.grid[, 2], func)
  grid.list <- list(
    x = xy.grid[, 1], y = xy.grid[, 2],
    z = z.grid
  )
  return(grid.list)
}

### plotPointProc
#' @noRd
".plotpointproc.Poisson" <- function(x, dev) {
  K <- x@K
  if (.check.grDevice() && dev) {
    dev.new(title = "Point Process Representation")
  }
  if (min(x@par$lambda) < 1) {
    lambda <- log(x@par$lambda)
  } else {
    lambda <- x@par$lambda
  }
  y.grid <- rep(0, K)
  size.grid <- as.vector(x@weight * 4)
  col.grid <- gray.colors(K,
    start = 0.2,
    end = 0.5
  )
  plot(lambda, y.grid,
    pch = 20, col = col.grid,
    cex = size.grid, cex.lab = .7, cex.axis = .7,
    main = "", ylab = "", xlab = ""
  )
  mtext(
    side = 1, bquote(lambda), cex = .7, cex.lab = .7,
    line = 3
  )
  legend.names <- list("", K)
  for (k in seq(1, K)) {
    legend.names[[k]] <- bquote(lambda[.(k)])
  }
  legend("topright",
    legend = do.call(expression, legend.names),
    col = col.grid, fill = col.grid
  )
}

### Has
### Checks if a 'model' object has specified parameters.
#' @noRd
".haspar.Model" <- function(obj, verbose) {
  if (length(obj@par) > 0) {
    dist <- obj@dist
    if (dist %in% c("poisson", "cond.poisson")) {
      .haspar.poisson.Model(obj, verbose)
    } else if (dist == "binomial") {
      .haspar.binomial.Model(obj, verbose)
    } else if (dist == "exponential") {
      .haspar.exponential.Model(obj, verbose)
    } else if (dist == "normal") {
      .haspar.normal.Model(obj, verbose)
    } else if (dist == "student") {
      .haspar.student.Model(obj, verbose)
    } else if (dist == "normult") {
      .haspar.normult.Model(obj, verbose)
    } else if (dist == "studmult") {
      .haspar.studmult.Model(obj, verbose)
    }
  } else {
    if (verbose) {
      stop(paste("Slot 'par' of 'model' object is ",
        "empty.",
        sep = ""
      ))
    } else {
      return(FALSE)
    }
  }
}

### -----------------------------------------------------------------
### .haspar.poisson.Mode
### @description    Checks if a Poisson model has fully specified
###                 parameters. If verbose is set to TRUE an error
###                 is thrown.
### @par    obj     an S4 object of class 'model'
### @par    verbose an object of class 'logical'
### @return         either TRUE or FALSE if parameters are fully
###                 specified or not. In case verbose == FALSE an
###                 error is thrown.
### -----------------------------------------------------------------
#' @noRd
".haspar.poisson.Model" <- function(obj, verbose) {
  if (length(obj@par) == 0) {
    if (verbose) {
      stop("Slot @par in 'model' object is empty.",
        call. = FALSE
      )
    } else {
      return(FALSE)
    }
  } else {
    if (!"lambda" %in% names(obj@par)) {
      if (verbose) {
        stop(paste("Wrong specification of slot @par ",
          "in 'model' object. Binomial models ",
          "need a parameter vector named 'lambda'.",
          sep = ""
        ), call. = FALSE)
      } else {
        return(FALSE)
      }
    } else {
      if (length(obj@par$lambda) != obj@K) {
        if (verbose) {
          stop(paste("Wrong specification of slot @par of ",
            "'model' object. Slot @K does not match ",
            "dimension of parameters in @par$lambda.",
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
### -------------------------------------------------------------------
### .haspar.binomial.Model
### @description    Checks if a Binomial model has fully specified
###                 parameters. If verbose is set to TRUE an error is
###                 thrown.
### @par    obj     an S4 object of class 'model'
### @par    verbose an object of class 'logical'
### @return         either TRUE or FALSE if parameters are fully
###                 specified or not. In case verbose == TRUE an
###                 error is thrown.
### -------------------------------------------------------------------
#' @noRd
".haspar.binomial.Model" <- function(obj, verbose) {
  if (length(obj@par) == 0) {
    if (verbose) {
      stop("Slot @par in 'model' object is empty.",
        call. = FALSE
      )
    } else {
      return(FALSE)
    }
  } else {
    if (!"p" %in% names(obj@par)) {
      if (verbose) {
        stop(paste("Wring specification of slot @par ",
          "in 'model' object. Binomial models ",
          "need a parameter named 'p'.",
          sep = ""
        ),
        call. = FALSE
        )
      } else {
        return(FALSE)
      }
    } else {
      if (length(obj@par$p) != obj@K) {
        if (verbose) {
          stop(paste("Wrong specification of slot @par of ",
            "'model' object. Slot @K does not ",
            "match the dimension of parameters ",
            "in @par$p.",
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

### ------------------------------------------------------------------
### .haspar.exponential.Model
### @description    Checks if an Exponential model has fully specified
###                 parameters. If verbose is set to TRUE an error is
###                 thrown.
### @param  obj     an S4 object of class 'model'
### @param  verbose an object of class 'logical'
### @return either TRUE or FALSE if parameters are fully specified or
###         nor. In case verbose == TRUE an error is thrown .
### ------------------------------------------------------------------
#' @noRd
".haspar.exponential.Model" <- function(obj, verbose) {
  if (length(obj@par) == 0) {
    if (verbose) {
      stop("Slot @par in 'model' object is empty",
        call. = FALSE
      )
    } else {
      return(FALSE)
    }
  } else {
    if (!"lambda" %in% names(obj@par)) {
      if (verbose) {
        stop(paste("Wrong specification of slot @par ",
          "in 'model' object. Exponential ",
          "models need a parameter named ",
          "'lambda'.",
          sep = ""
        ),
        call. = FALSE
        )
      } else {
        return(FALSE)
      }
    } else {
      if (length(obj@par$lambda) != obj@K) {
        if (verbose) {
          stop(paste("Wrong specification of slot @par in ",
            "'model' object. Number of Exponential ",
            "parameters in @par$lambda must match ",
            "number of components in slot @K.",
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

### ------------------------------------------------------------------
### .haspar.normal.Model
### @description    Checks if a Normal model has fully specified
###                 parameters. If verbose is set to TRUE an error is
###                 thrown.
### @param  obj     an S4 object of class 'model'
### @param  verbose an object of class 'logical'
### @return either TRUE or FALSE if parameters are fully specified or
###         not. In case verbose == TRUE an error is thrown .
### ------------------------------------------------------------------
#' @noRd
".haspar.normal.Model" <- function(obj, verbose) {
  K <- obj@K
  if (length(obj@par) == 0) {
    if (verbose) {
      stop("Slot @par in 'model' object is empty.",
        call. = FALSE
      )
    } else {
      return(FALSE)
    }
  } else {
    if (!("mu" %in% names(obj@par))) {
      if (verbose) {
        stop(paste("Wrong specification of slot @par ",
          "in 'model' object. Normal models ",
          "need a mean vector named 'mu'.",
          sep = ""
        ), call. = FALSE)
      } else {
        return(FALSE)
      }
    } else {
      if (length(obj@par$mu) != K) {
        if (verbose) {
          stop(paste("Wrong specification of slot @par ",
            "in 'model' object. Slot @K does ",
            "not match dimension of parameter ",
            "@par$mu.",
            sep = ""
          ), call. = FALSE)
        } else {
          return(FALSE)
        }
      } else {
        if (!("sigma" %in% names(obj@par))) {
          if (verbose) {
            stop(paste("Wrong specification of slot @par ",
              "in 'model' object. Normal models ",
              "need a standard deviation vector ",
              "named 'sigma'.",
              sep = ""
            ),
            call. = FALSE
            )
          } else {
            return(FALSE)
          }
        } else {
          if (length(obj@par$sigma) != K) {
            if (verbose) {
              stop(paste("Wrong specification of slot @par ",
                "in 'model' object. Slot @K does ",
                "not match dimension of parameter ",
                "par@$sigma.",
                sep = ""
              ), call. = FALSE)
            }
          } else {
            return(TRUE)
          }
        }
      }
    }
  }
}

### ------------------------------------------------------------------
### .haspar.normal.Model
### @description    Checks if a Normal model has fully specified
###                 parameters. If verbose is set to TRUE an error is
###                 thrown.
### @param  obj     an S4 object of class 'model'
### @param  verbose an object of class 'logical'
### @return either TRUE or FALSE if parameters are fully specified or
###         not. In case verbose == TRUE an error is thrown .
### ------------------------------------------------------------------
#' @noRd
".haspar.normult.Model" <- function(obj, verbose) {
  K <- obj@K
  if (length(obj@par) == 0) {
    if (verbose) {
      stop("Slot @par in 'model' object is empty.",
        call. = FALSE
      )
    } else {
      return(FALSE)
    }
  } else {
    if (!("mu" %in% names(obj@par))) {
      if (verbose) {
        stop(paste("Wrong specification of slot @par ",
          "in 'model' object. Normal models ",
          "need a mean vector named 'mu'.",
          sep = ""
        ), call. = FALSE)
      } else {
        return(FALSE)
      }
    } else {
      if (ncol(obj@par$mu) != K) {
        if (verbose) {
          stop(paste("Wrong specification of slot @par ",
            "in 'model' object. Slot @K does ",
            "not match dimension of parameter ",
            "@par$mu.",
            sep = ""
          ), call. = FALSE)
        } else {
          return(FALSE)
        }
      } else {
        if (!("sigma" %in% names(obj@par))) {
          if (verbose) {
            stop(paste("Wrong specification of slot @par ",
              "in 'model' object. Normal models ",
              "need a standard deviation vector ",
              "named 'sigma'.",
              sep = ""
            ),
            call. = FALSE
            )
          } else {
            return(FALSE)
          }
        } else {
          if (dim(obj@par$sigma)[3] != K) {
            if (verbose) {
              stop(paste("Wrong specification of slot @par ",
                "in 'model' object. Slot @K does ",
                "not match dimension of parameter ",
                "par@$sigma.",
                sep = ""
              ), call. = FALSE)
            }
          } else {
            return(TRUE)
          }
        }
      }
    }
  }
}

### ------------------------------------------------------------------
### .haspar.student.Model
### @description    Checks if a Normal model has fully specified
###                 parameters. If verbose is set to TRUE an error is
###                 thrown.
### @param  obj     an S4 object of class 'model'
### @param  verbose an object of class 'logical'
### @return either TRUE or FALSE if parameters are fully specified or
###         not. In case verbose == TRUE an error is thrown .
### ------------------------------------------------------------------
#' @noRd
".haspar.student.Model" <- function(obj, verbose) {
  K <- obj@K
  if (length(obj@par) == 0) {
    if (verbose) {
      stop("Slot @par in 'model' object is empty.",
        call. = FALSE
      )
    } else {
      return(FALSE)
    }
  } else {
    if (!("mu" %in% names(obj@par))) {
      if (verbose) {
        stop(paste("Wrong specification of slot @par ",
          "in 'model' object. Student-t models ",
          "need a mean vector named 'mu'.",
          sep = ""
        ), call. = FALSE)
      } else {
        return(FALSE)
      }
    } else {
      if (length(obj@par$mu) != K) {
        if (verbose) {
          stop(paste("Wrong specification of slot @par ",
            "in 'model' object. Slot @K does ",
            "not match dimension of parameter ",
            "@par$mu.",
            sep = ""
          ), call. = FALSE)
        } else {
          return(FALSE)
        }
      } else {
        if (!("sigma" %in% names(obj@par))) {
          if (verbose) {
            stop(paste("Wrong specification of slot @par ",
              "in 'model' object. Student-t models ",
              "need a standard deviation vector ",
              "named 'sigma'.",
              sep = ""
            ),
            call. = FALSE
            )
          } else {
            return(FALSE)
          }
        } else {
          if (length(obj@par$sigma) != K) {
            if (verbose) {
              stop(paste("Wrong specification of slot @par ",
                "in 'model' object. Slot @K does ",
                "not match dimension of parameter ",
                "par@$sigma.",
                sep = ""
              ), call. = FALSE)
            }
          } else {
            if (!"df" %in% names(obj@par)) {
              if (verbose) {
                stop(paste("Wrong specification of slot @par ",
                  "in 'model' object. Student-t models ",
                  "need a vector with degrees of freedom ",
                  "named 'df'.",
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

### ------------------------------------------------------------------
### .haspar.student.Model
### @description    Checks if a Normal model has fully specified
###                 parameters. If verbose is set to TRUE an error is
###                 thrown.
### @param  obj     an S4 object of class 'model'
### @param  verbose an object of class 'logical'
### @return either TRUE or FALSE if parameters are fully specified or
###         not. In case verbose == TRUE an error is thrown .
### ------------------------------------------------------------------
#' @noRd
".haspar.studmult.Model" <- function(obj, verbose) {
  K <- obj@K
  if (length(obj@par) == 0) {
    if (verbose) {
      stop("Slot @par in 'model' object is empty.",
        call. = FALSE
      )
    } else {
      return(FALSE)
    }
  } else {
    if (!("mu" %in% names(obj@par))) {
      if (verbose) {
        stop(paste("Wrong specification of slot @par ",
          "in 'model' object. Student-t models ",
          "need a mean vector named 'mu'.",
          sep = ""
        ), call. = FALSE)
      } else {
        return(FALSE)
      }
    } else {
      if (ncol(obj@par$mu) != K) {
        if (verbose) {
          stop(paste("Wrong specification of slot @par ",
            "in 'model' object. Slot @K does ",
            "not match dimension of parameter ",
            "@par$mu.",
            sep = ""
          ), call. = FALSE)
        } else {
          return(FALSE)
        }
      } else {
        if (!("sigma" %in% names(obj@par))) {
          if (verbose) {
            stop(paste("Wrong specification of slot @par ",
              "in 'model' object. Student-t models ",
              "need a standard deviation vector ",
              "named 'sigma'.",
              sep = ""
            ),
            call. = FALSE
            )
          } else {
            return(FALSE)
          }
        } else {
          if (dim(obj@par$sigma)[3] != K) {
            if (verbose) {
              stop(paste("Wrong specification of slot @par ",
                "in 'model' object. Slot @K does ",
                "not match dimension of parameter ",
                "par@$sigma.",
                sep = ""
              ), call. = FALSE)
            }
          } else {
            if (!"df" %in% names(obj@par)) {
              if (verbose) {
                stop(paste("Wrong specification of slot @par ",
                  "in 'model' object. Student-t models ",
                  "need a vector with degrees of freedom ",
                  "named 'df'.",
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

### Validity
### Validity checking of model objects is implemented
### in two versions: an initializing version relying partly
### on warnings and amore restrictive version relying exclusively
### on errors.
### The less restrictive validity check is used in setters and
### and the fully restrictive version in the constructor and later
### usage of model object (e.g. see 'mcmcstart()')
### -----------------------------------------------------------------------------
### .init.valid.Model
### @description    Initial validity check for model object
### @par    obj     a model object
### @return         An error in case certain conditions are failed or there are
###                 inconsistencies.
### @see            ?model, ?vignette('finmix'), .init.valid.*, .valid.*
### @author         Lars Simon Zehnder
### -----------------------------------------------------------------------------
#' @noRd
".init.valid.Model" <- function(obj) {
  .valid.dist.Model(obj)
  .init.valid.K.Model(obj)
  .init.valid.r.Model(obj)
  .init.valid.par.Model(obj)
  .init.valid.weight.Model(obj)
  .init.valid.T.Model(obj)
}

### -----------------------------------------------------------------------------
### .init.Model
### @description    Validity check for model object
### @par    obj     a model object
### @return         An error in case certain conditions are failed or a warning
###                 if there are inconsistencies.
### @see            ?model, ?vignette('finmix'), .init.valid.*, .valid.*
### @author         Lars Simon Zehnder
### -----------------------------------------------------------------------------
#' @noRd
".valid.Model" <- function(obj) {
  .valid.dist.Model(obj)
  .valid.K.Model(obj)
  .valid.r.Model(obj)
  .valid.par.Model(obj)
  .valid.weight.Model(obj)
  .valid.T.Model(obj)
}

### ----------------------------------------------------------------------------
### .valid.dist.Model
### @description    Initial validity check for the distribution of a finite
###                 mixture model
### @par    obj     a model object
### @return         An error in case the distribution is unknown.
### @see            ?model, ?vignette('finmix')i
### ----------------------------------------------------------------------------
#' @noRd
".valid.dist.Model" <- function(obj) {
  dists <- c(
    "normal", "normult", "exponential",
    "student", "studmult", "poisson",
    "cond.poisson", "binomial"
  )
  indicmod.dists <- c("multinomial")
  if (length(obj@dist) > 0) {
    if (!(obj@dist %in% dists)) {
      stop(paste("Unknown distribution in slot 'dist' ",
        "of 'model' object.",
        sep = ""
      ),
      call. = FALSE
      )
    } else {
      if (!(obj@indicmod %in% indicmod.dists)) {
        stop(paste("Unknown indicator distribution in slot ",
          "'indicmod' of 'model' object.",
          sep = ""
        ),
        call. = FALSE
        )
      }
    }
  }
}

### ----------------------------------------------------------------------------
### .init.valid.K.Model
### @description    Initial validity check for the number of components K of
###                 a finite mixture model.
### @par    obj     a model object
### @return         An error if the number of components are not a positive
###                 integer
### @see            ?model, ?vignette('finmix')
### @author         Lars Simon Zehnder
### ----------------------------------------------------------------------------
#' @noRd
".init.valid.K.Model" <- function(obj) {
  if (obj@K < 1) {
    stop(paste("Wrong specification of slot 'K' of ",
      "'model' object. Number of components ",
      "must be a positive integer.",
      sep = ""
    ),
    call. = FALSE
    )
  } else {
    if (!all(is.na(obj@weight))) {
      if (obj@K != ncol(obj@weight)) {
        stop(paste("Dimension of slot 'weight' in ",
          "'model' object does not match ",
          "number of components in slot 'K'.",
          sep = ""
        ),
        call. = FALSE
        )
      }
    }
    .init.valid.par.Model(obj)
  }
}

### ----------------------------------------------------------------------------
### .valid.K.Model
### @description    Validity check for the number of components K of
###                 a finite mixture model.
### @par    obj     a model object
### @return         An error if the number of components are not a positive
###                 integer and a warning if the number of components do not
###                 match the dimension of the weights.
### @see            ?model, ?vignette('finmix')
### @author         Lars Simon Zehnder
### ----------------------------------------------------------------------------
#' @noRd
".valid.K.Model" <- function(obj) {
  if (obj@K < 1) {
    stop(paste("Wrong specification of slot 'K' of ",
      "'model' object. Number of components ",
      "must be a positive integer.",
      sep = ""
    ),
    call. = FALSE
    )
  } else {
    if (!all(is.na(obj@weight))) {
      if (obj@K != ncol(obj@weight)) {
        warning(paste("Dimension of slot 'weight' in ",
          "'model' object does not match ",
          "number of components in slot 'K'.",
          sep = ""
        ),
        call. = FALSE
        )
      }
    }
    .valid.par.Model(obj)
  }
}

### ----------------------------------------------------------------------------
### .init.valid.r.Model
### @description    Initial validity check for variable dimension r.
### @par    obj     a model object
### @return         An error in case the variable dimension r is not a positive
###                 integer or the dimension does not fit the distribution model.
### @see            ?model, ?vignette('finmix')
### @author         Lars Simon Zehnder
### ----------------------------------------------------------------------------
#' @noRd
".init.valid.r.Model" <- function(obj) {
  univ <- .get.univ.Model()
  multiv <- .get.multiv.Model()
  if (obj@r < 1) {
    stop(paste("Wrong specification of slot 'r' ",
      "in 'model' object. Dimension of ",
      "variables must be a positive integer.",
      sep = ""
    ),
    call. = FALSE
    )
  } else {
    if ((obj@dist %in% univ) && obj@r > 1) {
      stop(paste("Wrong specification of slot 'r' ",
        "in 'model' object. Univariate ",
        "distributions can only have one ",
        "dimension.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if ((obj@dist %in% multiv) && obj@r < 2) {
      stop(paste("Wrong specification of slot 'r' ",
        "in 'model' object. Multivariate ",
        "distributions must have dimension ",
        "greater one.",
        sep = ""
      ),
      call. = FALSE
      )
    }
  }
}

### ----------------------------------------------------------------------------
### .init.valid.r.Model
### @description    Initial validity check for variable dimension r.
### @par    obj     a model object
### @return         An error in case the variable dimension r is not a positive
###                 integer or a warning if the dimension does not fit the
###                 distribution model.
### @see            ?model, ?vignette('finmix')
### @author         Lars Simon Zehnder
### ----------------------------------------------------------------------------
#' @noRd
".valid.r.Model" <- function(obj) {
  univ <- .get.univ.Model()
  multiv <- .get.multiv.Model()
  if (obj@r < 1) {
    stop(paste("Wrong specification of slot 'r' ",
      "in 'model' object. Dimension of ",
      "variables must be positive.",
      sep = ""
    ),
    call. = FALSE
    )
  } else {
    if ((obj@dist %in% univ) && obj@r > 1) {
      stop(paste("Wrong specification of slot 'r' ",
        "in 'model' object. Univariate ",
        "distributions can only have one ",
        "dimension.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if ((obj@dist %in% multiv) && obj@r < 2) {
      stop(paste("Wrong specification of slot 'r' ",
        "in 'model' object. Multivariate ",
        "distributions must have dimension ",
        "greater one.",
        sep = ""
      ),
      call. = FALSE
      )
    }
  }
}

### ----------------------------------------------------------------------------
### .init.valid.weight.Model
### @description    Initial validity check for the weights of a finite mixture
###                 model.
### @par    obj     a model object
### @return         An error if the dimension of the weight vector does not fit
###                 the model or if the weights do not sum to 1, are negative or
###                 larger than one.
### @see            ?model, ?vignette('finmix')
### @author         Lars Simon Zehnder
### ----------------------------------------------------------------------------
#' @noRd
".init.valid.weight.Model" <- function(obj) {
  if (!all(is.na(obj@weight))) {
    if (nrow(obj@weight) > 1) {
      stop(paste("Wrong dimension of slot 'weight' in ",
        "'model' object. Dimension of slot ",
        "'weight' must be 1 x K.",
        sep = ""
      ),
      call. = FALSE
      )
    } else {
      if (ncol(obj@weight) != obj@K) {
        stop(paste("Wrong number of weights in slot 'weight' of ",
          "'model' object. Number of weights does not ",
          "match number of components in slot 'K'.",
          sep = ""
        ),
        call. = FALSE
        )
      } else {
        if (is.integer(obj@weight)) {
          stop(paste("Wrong specification of slot 'weight' of ",
            "'model' object. Weights must be of type ",
            "'numeric'.",
            sep = ""
          ),
          call. = FALSE
          )
        }
        if (!is.numeric(obj@weight)) {
          stop(paste("Wrong specification of slot 'weight' of ",
            "'model' object. Weights must be of type ",
            "'numeric'.",
            sep = ""
          ),
          call. = FALSE
          )
        }
        if (any(obj@weight <= 0) || any(obj@weight >= 1)) {
          stop(paste("Weights in slot 'weight' of 'model' ",
            "object must be positive.",
            sep = ""
          ),
          call. = FALSE
          )
        } else {
          if (round(sum(obj@weight)) != 1) {
            stop(paste("Weights in slot 'weight' of 'model' ",
              "object must sum to one.",
              sep = ""
            ),
            call. = FALSE
            )
          }
        }
      }
    }
  }
}

### ------------------------------------------------------------------------------------
### .valid.weight.Model
### @description    Validity check for the weights of a finite mixture model.
### @par    obj     a model object
### @return         An error if the weights are not of type 'numeric' and a warning
###                 if the weigths do not conform to the number of components K,
###                 do not sum to one or are not values between 0 and 1.
### @see            ?model, ?vignette('finmix')
### @author         Lars Simon Zehnder
### -------------------------------------------------------------------------------------
#' @noRd
".valid.weight.Model" <- function(obj) {
  if (!all(is.na(obj@weight))) {
    if (nrow(obj@weight) > 1) {
      warning(paste("Wrong dimension of slot 'weight' in ",
        "'model' object. Dimension of slot ",
        "'weight' must be 1 x K.",
        sep = ""
      ),
      call. = FALSE
      )
    } else {
      if (ncol(obj@weight) != obj@K) {
        warning(paste("Wrong number of weights in slot 'weight' of ",
          "'model' object. Number of weights does not ",
          "match number of components in slot 'K'.",
          sep = ""
        ),
        call. = FALSE
        )
      } else {
        if (is.integer(obj@weight)) {
          stop(paste("Wrong specification of slot 'weight' of ",
            "'model' object. Weights must be of type ",
            "'numeric'.",
            sep = ""
          ),
          call. = FALSE
          )
        }
        if (!is.numeric(obj@weight)) {
          stop(paste("Wrong specification of slot 'weight' of ",
            "'model' object. Weights must be of type ",
            "'numeric'.",
            sep = ""
          ),
          call. = FALSE
          )
        }
        if (any(obj@weight <= 0) || any(obj@weight >= 1)) {
          warning(paste("Weights in slot 'weight' of 'model' ",
            "object must be positive.",
            sep = ""
          ),
          call. = FALSE
          )
        } else {
          if (round(sum(obj@weight)) != 1) {
            warning(paste("Weights in slot 'weight' of 'model' ",
              "object must sum to one.",
              sep = ""
            ),
            call. = FALSE
            )
          }
        }
      }
    }
  }
}

### -------------------------------------------------------------------------------------
### .init.valid.T.Model
### @description    Initial validity check for the repetitions of a Binomial mixture.
### @par    obj     a model object
### @return         An error in case the reptitions are not of type integer, have
###                 the wrong dimension, or non-positive values.
### @see            ?model, ?vignette('finmix')
### --------------------------------------------------------------------------------------
#' @noRd
".init.valid.T.Model" <- function(obj) {
  if (!all(is.na(obj@T))) {
    if (!is.integer(obj@T)) {
      stop(paste("Wrong type of slot 'T' in 'model' object ",
        "Repetitions must be of type 'integer'.",
        sep = ""
      ),
      call. = FALSE
      )
    }
    if (nrow(obj@T) > 1 && ncol(obj@T) > 1) {
      stop(paste("Wrong dimension of slot 'T' in 'model' ",
        "object. Repetitions can only be ",
        "one-dimensional",
        sep = ""
      ),
      call. = FALSE
      )
    }
    if (any(obj@T < 1)) {
      stop(paste("Wrong specification of slot 'T' in 'model' ",
        "object. Repetitions must be positive integers ",
        "or NA.",
        sep = ""
      ),
      call. = FALSE
      )
    }
  }
}

### -------------------------------------------------------------------------------------
### .valid.T.Model
### @description    Validity check for the repetitions of a Binomial mixture.
### @par    obj     a model object
### @return         An error in case the reptitions are not of type integer, have
###                 the wrong dimension, or non-positive values.
### @see            ?model, ?vignette('finmix')
### --------------------------------------------------------------------------------------
#' @noRd
".valid.T.Model" <- function(obj) {
  if (!all(is.na(obj@T))) {
    if (!is.integer(obj@T)) {
      stop(paste("Wrong type of slot 'T' in 'model' object ",
        "Repetitions must be of type 'integer'.",
        sep = ""
      ),
      call. = FALSE
      )
    }
    if (nrow(obj@T) > 1 && ncol(obj@T) > 1) {
      stop(paste("Wrong dimension of slot 'T' in 'model' ",
        "object. Repetitions can only be ",
        "one-dimensional",
        sep = ""
      ),
      call. = FALSE
      )
    }
    if (any(obj@T < 1)) {
      stop(paste("Wrong specification of slot 'T' in 'model' ",
        "object. Repetitions must be positive integers ",
        "or NA.",
        sep = ""
      ),
      call. = FALSE
      )
    }
  }
}


### -------------------------------------------------------------------------------
### .init.valid.par.Model
### @description    Initial validity check of model parameters
### @par    obj     a model object
### @return         An error if parameters fail certain conditions
### @detail         This validity check is called in the S4 constructor
###                 'model()' and ensures that the user constructs an inherently
###                 consistent model object.
### @see            ?model, ?vignette('finmix')
### @author         Lars Simon Zehnder
### --------------------------------------------------------------------------------
#' @noRd
".init.valid.par.Model" <- function(obj) {
  dist <- obj@dist
  if (length(obj@par) > 0) {
    if (dist %in% c("poisson", "cond.poisson")) {
      .init.valid.Poisson.Model(obj)
    } else if (dist == "binomial") {
      .init.valid.Binomial.Model(obj)
    } else if (dist == "normal") {
      .init.valid.Normal.Model(obj)
    } else if (dist == "normult") {
      .init.valid.Normult.Model(obj)
    } else if (dist == "student") {
      .init.valid.Student.Model(obj)
    } else if (dist == "studmult") {
      .init.valid.Studmult.Model(obj)
    }
  }
}

### -------------------------------------------------------------------------------
### .valid.par.Model
### @description    Validity check of model parameters
### @par    obj     a model object
### @return         An error if parameters fail certain necessary conditions and
###                 a warning if parameters fail consistency.
### @detail         This validity check is called in the setters to ensure that
###                 slots can be changed without errors but help the user to
###                 end up with an inherently consistent model object.
### @see            ?model, ?vignette('finmix')
### @author         Lars Simon Zehnder
### --------------------------------------------------------------------------------
#' @noRd
".valid.par.Model" <- function(obj) {
  dist <- obj@dist
  if (length(obj@par) > 0) {
    if (dist %in% c("poisson", "cond.poisson")) {
      .valid.Poisson.Model(obj)
    } else if (dist == "binomial") {
      .valid.Binomial.Model(obj)
    } else if (dist == "exponential") {
      .valid.Exponential.Model(obj)
    } else if (dist == "normal") {
      .valid.Normal.Model(obj)
    } else if (dist == "normult") {
      .valid.Normult.Model(obj)
    } else if (dist == "student") {
      .valid.Student.Model(obj)
    } else if (dist == "studmult") {
      .valid.Studmult.Model(obj)
    }
  }
}

### -----------------------------------------------------------------------------
### .init.valid.Poisson.Model
### @description    Initial validity check for parameters of a Poisson mixture.
### @par    obj     a model object
### @return         An error if parameters fail certain conditions.
### @detail     This initial validity check is called in the S4 constructor
###             'model()' and ensures that the user constructs an inherently
###             consistent model object.
###             The parameter list must contain an element 'lambda' that is
###             an 1 x K array, vector or matrix with numeric or integer values
###             all positive.
### @see        ?model
### @author     Lars Simon Zehnder
### -------------------------------------------------------------------------------
#' @noRd
".init.valid.Poisson.Model" <- function(obj) {
  if (length(obj@par) > 0) {
    if ("lambda" %in% names(obj@par)) {
      if (!is.array(obj@par$lambda) && !is.vector(obj@par$lambda) &&
        !is.matrix(obj@par$lambda)) {
        stop(paste("Wrong specification of slot @par: ",
          "Poisson parameters must be either an ",
          "array, a vector or a matrix of dimension ",
          "1 x K.",
          sep = ""
        ),
        call. = FALSE
        )
      }
      obj@par$lambda <- as.vector(obj@par$lambda)
      if (!is.numeric(obj@par$lambda) && !is.integer(obj@par$lambda)) {
        stop(paste("Wrong specification in slot 'par' of 'model' object. ",
          "Parameters must be of type 'numeric' or 'integer'.",
          sep = ""
        ),
        call. = FALSE
        )
      }
      if (length(obj@par$lambda) != obj@K) {
        warning(paste("Wrong specification of slot @par: ",
          "lambda must be either an array, a vector ",
          "or a matrix of dimension 1 x K.",
          sep = ""
        ),
        call. = FALSE
        )
      } else {
        if (any(obj@par$lambda <= 0)) {
          stop(paste("Wrong specification of slot @par: ",
            "Poisson parameters ",
            "must be positive.",
            sep = ""
          ),
          call. = FALSE
          )
        }
      }
    } else {
      warning(paste("Wrong specification of slot 'par' in 'model' object. ",
        "Poisson parameters must be named 'lambda'.",
        sep = ""
      ),
      call. = FALSE
      )
    }
  }
}

### -----------------------------------------------------------------------------
### .valid.Poisson.Model
### @description    Validity check for parameters of a Poisson mixture.
### @par    obj     a model object
### @return         An error if parameters do fail certain necessary conditions.
###                 A warning if parameters do fail consistency.
### @detail     This validity check is called in the setters to ensure that
###             slots can be changed without errors but help the user to
###             get a inherently consistent model object.
###             The parameter list must contain an element 'lambda' that is
###             an 1 x K array, vector or matrix with numeric or integer values
###             all positive.
### @see        $model
### @author     Lars Simon Zehnder
### -----------------------------------------------------------------------------
#' @noRd
".valid.Poisson.Model" <- function(obj) {
  if (length(par) > 0) {
    if ("lambda" %in% names(obj@par)) {
      if (!is.array(obj@par$lambda) && !is.vector(obj@par$lambda) &&
        !is.matrix(obj@par$lambda)) {
        stop(paste("Wrong specification of slot @par: ",
          "Poisson parameters must be either an ",
          "array, a vector or a matrix of dimension ",
          "1 x K.",
          sep = ""
        ),
        call. = FALSE
        )
      }
      obj@par$lambda <- as.vector(obj@par$lambda)
      if (!is.numeric(obj@par$lambda) && !is.integer(obj@par$lambda)) {
        stop(paste("Wrong specification in slot 'par' of 'model' object. ",
          "Parameters must be of type 'numeric' or 'integer'.",
          sep = ""
        ),
        call. = FALSE
        )
      }
      if (length(obj@par$lambda) != obj@K) {
        warning(paste("Wrong specification of slot @par: ",
          "lambda must be either an array, a vector ",
          "or a matrix of dimension 1 x K.",
          sep = ""
        ),
        call. = FALSE
        )
      } else {
        if (any(obj@par$lambda <= 0)) {
          stop(paste("Wrong specification of slot @par: ",
            "Poisson parameters ",
            "must be positive.",
            sep = ""
          ),
          call. = FALSE
          )
        }
      }
    } else {
      stop(paste("Wrong specification of slot 'par' in 'model' object. ",
        "Poisson parameters must be named 'lambda'.",
        sep = ""
      ),
      call. = FALSE
      )
    }
  }
}

### ------------------------------------------------------------------------------
### .init.valid.Binomial.Model
### @description    Initial validity check for parameters of a Binomial mixture.
### @par    obj     a model object
### @return         An error if parameters fail certain conditions
### @detail         This initial validity check is called in the S4 constructor
###                 'model()' and ensures that the user constructs an inherently
###                 consistent model object.
###                 The parameter list must contain an 1 x K array, vector, or
###                 matrix with probabilities, all between 0 and 1.
### @see            ?model, ?vignette('finmix')
### @author         Lars Simon Zehnder
### ------------------------------------------------------------------------------
#' @noRd
".init.valid.Binomial.Model" <- function(obj) {
  if (length(obj@par)) {
    if (!"p" %in% names(obj@par)) {
      stop(paste("Wrong specification of slot @par: ",
        "Binomial mixtures need a ",
        "probability vector named 'p'.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!is.array(obj@par$p) && !is.vector(obj@par$p) &&
      !is.matrix(obj@par$p)) {
      stop(paste("Wrong specification of slot @par: ",
        "p must be either an array, a vector ",
        "or a matrix of dimension 1 x K",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!all(is.numeric(obj@par$p) || is.integer(obj@par$p))) {
      stop(paste("Wrong specification of slot @par: ",
        "parameters must be either of type ,",
        "'numeric' or 'integer'.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (length(obj@par$p) != obj@K) {
      stop(paste("Wrong specification of slot @par: ",
        "p must be an array, a vector ",
        "or a matrix of dimension 1 x K",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!all(obj@par$p > 0 && obj@par$p < 1)) {
      stop(paste("Wrong specification of slot @par: ",
        "Binomial parameters must be all ",
        "between 0 and 1.",
        sep = ""
      ),
      call. = FALSE
      )
    }
  }
  if (dim(obj@T)[1] > 1 && dim(obj@T)[2] > 1) {
    stop(paste(
      "Dimensions of repetitions 'T' for binomial mixture",
      "model do not match conditions. Only one-dimensional",
      "repetitions can be used in a binomial mixture model."
    ), sep = "")
  }
}

### ------------------------------------------------------------------------------
### .valid.Binomial.Model
### @description    Validity check for parameters of a Binomial mixture.
### @par    obj     a model object
### @return         An error if parameters fail certain necessary conditions and
###                 a warning if parameters fail consistency
### @detail         This validity check is called in the setters to ensure that
## ä                 slots can be changed without errors but help the user to
###                 end up with an inherently consistent model object.
###                 The parameter list must contain an 1 x K array, vector, or
###                 matrix with probabilities, all between 0 and 1.
### @see            ?model, ?vignette('finmix')
### @author         Lars Simon Zehnder
### ------------------------------------------------------------------------------
#' @noRd
".valid.Binomial.Model" <- function(obj) {
  if (length(obj@par)) {
    if (!"p" %in% names(obj@par)) {
      warning(paste("Wrong specification of slot @par: ",
        "Binomial mixtures need a ",
        "probability vector named 'p'.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!is.array(obj@par$p) && !is.vector(obj@par$p) &&
      !is.matrix(obj@par$p)) {
      stop(paste("Wrong specification of slot @par: ",
        "p must be either an array, a vector ",
        "or a matrix of dimension 1 x K",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!all(is.numeric(obj@par$p) || is.integer(obj@par$p))) {
      stop(paste("Wrong specification of slot @par: ",
        "parameters must be either of type ,",
        "'numeric' or 'integer'.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (length(obj@par$p) != obj@K) {
      warning(paste("Wrong specification of slot @par: ",
        "p must be an array, a vector ",
        "or a matrix of dimension 1 x K",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!all(obj@par$p > 0 && obj@par$p < 1)) {
      stop(paste("Wrong specification of slot @par: ",
        "Binomial parameters must be all ",
        "between 0 and 1.",
        sep = ""
      ),
      call. = FALSE
      )
    }
  }
  if (dim(obj@T)[1] > 1 && dim(obj@T)[2] > 1) {
    stop(paste(
      "Dimensions of repetitions 'T' for binomial mixture",
      "model do not match conditions. Only one-dimensional",
      "repetitions can be used in a binomial mixture model."
    ), sep = "")
  }
}


### -----------------------------------------------------------------------------
### .init.valid.Exponential.Model
### @description    Initial validity check for parameters of a Exponential
###                 mixture.
### @par    obj     a model object
### @return         An error if parameters fail certain conditions.
### @detail     This initial validity check is called in the S4 constructor
###             'model()' and ensures that the user constructs an inherently
###             consistent model object.
###             The parameter list must contain an element 'lambda' that is
###             an 1 x K array, vector or matrix with numeric or integer values
###             all positive.
### @see        ?model
### @author     Lars Simon Zehnder
### -------------------------------------------------------------------------------
#' @noRd
".init.valid.Exponential.Model" <- function(obj) {
  if (length(obj@par) > 0) {
    if ("lambda" %in% names(obj@par)) {
      if (!is.array(obj@par$lambda) && !is.vector(obj@par$lambda) &&
        !is.matrix(obj@par$lambda)) {
        stop(paste("Wrong specification of slot @par: ",
          "Exponential parameters must be either an ",
          "array, a vector or a matrix of dimension ",
          "1 x K.",
          sep = ""
        ),
        call. = FALSE
        )
      }
      obj@par$lambda <- as.vector(obj@par$lambda)
      if (!is.numeric(obj@par$lambda) && !is.integer(obj@par$lambda)) {
        stop(paste("Wrong specification in slot 'par' of 'model' object. ",
          "Parameters must be of type 'numeric' or 'integer'.",
          sep = ""
        ),
        call. = FALSE
        )
      }
      if (length(obj@par$lambda) != obj@K) {
        warning(paste("Wrong specification of slot @par: ",
          "lambda must be either an array, a vector ",
          "or a matrix of dimension 1 x K.",
          sep = ""
        ),
        call. = FALSE
        )
      } else {
        if (any(obj@par$lambda <= 0)) {
          stop(paste("Wrong specification of slot @par: ",
            "Exponential parameters ",
            "must be positive.",
            sep = ""
          ),
          call. = FALSE
          )
        }
      }
    } else {
      warning(paste("Wrong specification of slot 'par' in 'model' object. ",
        "Exponential parameters must be named 'lambda'.",
        sep = ""
      ),
      call. = FALSE
      )
    }
  }
}

### -----------------------------------------------------------------------------
### .valid.Exponential.Model
### @description    Validity check for parameters of a Exponential mixture.
### @par    obj     a model object
### @return         An error if parameters do fail certain necessary conditions.
###                 A warning if parameters do fail consistency.
### @detail     This validity check is called in the setters to ensure that
###             slots can be changed without errors but help the user to
###             get a inherently consistent model object.
###             The parameter list must contain an element 'lambda' that is
###             an 1 x K array, vector or matrix with numeric or integer values
###             all positive.
### @see        $model
### @author     Lars Simon Zehnder
### -----------------------------------------------------------------------------
#' @noRd
".valid.Exponential.Model" <- function(obj) {
  if (length(par) > 0) {
    if ("lambda" %in% names(obj@par)) {
      if (!is.array(obj@par$lambda) && !is.vector(obj@par$lambda) &&
        !is.matrix(obj@par$lambda)) {
        stop(paste("Wrong specification of slot @par: ",
          "Exponential parameters must be either an ",
          "array, a vector or a matrix of dimension ",
          "1 x K.",
          sep = ""
        ),
        call. = FALSE
        )
      }
      obj@par$lambda <- as.vector(obj@par$lambda)
      if (!is.numeric(obj@par$lambda) && !is.integer(obj@par$lambda)) {
        stop(paste("Wrong specification in slot 'par' of 'model' object. ",
          "parameters must be of type 'numeric' or 'integer'.",
          sep = ""
        ),
        call. = FALSE
        )
      }
      if (length(obj@par$lambda) != obj@K) {
        warning(paste("Wrong specification of slot @par: ",
          "lambda must be either an array, a vector ",
          "or a matrix of dimension 1 x K.",
          sep = ""
        ),
        call. = FALSE
        )
      } else {
        if (any(obj@par$lambda <= 0)) {
          stop(paste("Wrong specification of slot @par: ",
            "Exponential parameters ",
            "must be positive.",
            sep = ""
          ),
          call. = FALSE
          )
        }
      }
    } else {
      stop(paste("Wrong specification of slot 'par' in 'model' object. ",
        "Exponential parameters must be named 'lambda'.",
        sep = ""
      ),
      call. = FALSE
      )
    }
  }
}

### ------------------------------------------------------------------------------
### .init.valid.Normal.Model
### @description    Initial validity check for parameters of a univariate
###                 Normal mixture.
### @par    obj     a model object
### @return         An error if parameters fail certain conditions
### @detail         This initial validity check is called in the S4 constructor
###                 'model()' and ensures that the user constructs an inherently
###                 consistent model object.
###                 The parameter list must contain the following elements:
###                     mu:     an 1 x K array, vector or matrix containing
###                             'numeric' or 'integer' values
###                     sigma:  an 1 x K array, vector or matrix containing
###                             'numeric' or 'integer' values, all positive
###                     df:     an 1 x K array, vector or matrix containing
###                             'numeric' or 'integer' values, all positive
### @see            ?model, ?vignette('finmix')
### @author         Lars Simon Zehnder
### -------------------------------------------------------------------------------
#' @noRd
".init.valid.Normal.Model" <- function(obj) {
  if (length(obj@par) > 0) {
    if (!"mu" %in% names(obj@par)) {
      stop(paste("Wrong specification of slot @par: ",
        "univariate Normal mixtures need ",
        "a mean matrix named 'mu'.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!is.array(obj@par$mu) && !is.vector(obj@par$mu) &&
      !is.matrix(obj@par$mu)) {
      stop(paste("Wrong specification of slot @par: ",
        "mu must be either an array, a vector ",
        "or a matrix of dimension 1 x K. ",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!all(is.numeric(as.vector(obj@par$mu)) ||
      is.integer(as.vector(obj@par$mu)))) {
      stop(paste("Wrong specification of slot @par: ",
        "Parameters must be of type 'numeric' or ",
        "'integer'.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (length(obj@par$mu) != obj@K) {
      stop(paste("Wrong specification of slot @par: ",
        "mu must be a matrix of dimension 1 x K ",
        "or a vector of size K.",
        sep = ""
      ),
      call. = FALSE
      )
    }
    if (!"sigma" %in% names(obj@par)) {
      stop(paste("Wrong specification of slot @par: ",
        "univariate Normal mictures need ",
        "a variance vector named ",
        "'sigma'",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!all(is.numeric(as.vector(obj@par$sigma)) ||
      is.integer(as.vector(obj@par$sigma)))) {
      stop(paste("Wrong specification of slot @par: ",
        "Parameters must be of type 'numeric' or ",
        "'integer'.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (any(obj@par$sigma <= 0)) {
      stop(paste("Wrong specification of slot @par: ",
        "sigma must contain variances, all ",
        "positive.",
        sep = ""
      ),
      .call = FALSE
      )
    } else if (!is.array(obj@par$sigma) && !is.vector(obj@par$sigma) &&
      !is.matrix(obj@par$sigma)) {
      stop(paste("Wrong specification of slot @par: ",
        "sigma must be either an array, a vector, ",
        "or a matrix of dimension 1 x K.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (length(obj@par$sigma) != obj@K) {
      stop(paste("Wrong specification of slot @par: ",
        "sigma must be either an array, a vector, ",
        "or a matrix, ",
        "or a matrix of dimension ",
        "1 x K.",
        sep = ""
      ),
      call. = FALSE
      )
    }
  }
}

### ------------------------------------------------------------------------------
### .valid.Normal.Model
### @description    Validity check for parameters of a univariate Normal
###                 mixture.
### @par    obj     a model object
### @return         An error if parameters fail certain necessary conditions and
###                 a warning if parameters fail consistency.
### @detail         This validity check is called in the setters to ensure that
###                 slots can be changed without errors but help the user to
###                 end up with an inherently consistent model object.
###                 The parameter list must contain the following elements:
###                     mu:     an 1 x K array, vector or matrix containing
###                             'numeric' or 'integer' values
###                     sigma:  an 1 x K array, vector or matrix containing
###                             'numeric' or 'integer' values, all positive
### @see            ?model, ?vignette('finmix')
### @author         Lars Simon Zehnder
### -------------------------------------------------------------------------------
#' @noRd
".valid.Normal.Model" <- function(obj) {
  if (length(obj@par) > 0) {
    if (!"mu" %in% names(obj@par)) {
      warning(paste("Wrong specification of slot @par: ",
        "univariate Normal mixtures need ",
        "a mean matrix named 'mu'.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!is.array(obj@par$mu) && !is.vector(obj@par$mu) &&
      !is.matrix(obj@par$mu)) {
      warning(paste("Wrong specification of slot @par: ",
        "mu must be either an array, a vector ",
        "or a matrix of dimension 1 x K. ",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!all(is.numeric(as.vector(obj@par$mu)) ||
      is.integer(as.vector(obj@par$mu)))) {
      stop(paste("Wrong specification of slot @par: ",
        "Parameters must be of type 'numeric' or ",
        "'integer'.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (length(obj@par$mu) != obj@K) {
      warning(paste("Wrong specification of slot @par: ",
        "mu must be a matrix of dimension 1 x K ",
        "or a vector of size K.",
        sep = ""
      ),
      call. = FALSE
      )
    }
    if (!"sigma" %in% names(obj@par)) {
      warning(paste("Wrong specification of slot @par: ",
        "univariate Normal mixtures need ",
        "a variance vector named ",
        "'sigma'",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!all(is.numeric(as.vector(obj@par$sigma)) ||
      is.integer(as.vector(obj@par$sigma)))) {
      stop(paste("Wrong specification of slot @par: ",
        "Parameters must be of type 'numeric' or ",
        "'integer'.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (any(obj@par$sigma <= 0)) {
      stop(paste("Wrong specification of slot @par: ",
        "sigma must contain variances, all ",
        "positive.",
        sep = ""
      ),
      .call = FALSE
      )
    } else if (!is.array(obj@par$sigma) && !is.matrix(obj@par$sigma) &&
      !is.matrix(obj@par$sigma)) {
      warning(paste("Wrong specification of slot @par: ",
        "sigma must be either an array, a vector, ",
        "or a matrix of dimension 1 x K.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (length(obj@par$sigma) != obj@K) {
      warning(paste("Wrong specification of slot @par: ",
        "sigma must be either an array, a vector, ",
        "or a matrix, ",
        "or a matrix of dimension ",
        "1 x K.",
        sep = ""
      ),
      call. = FALSE
      )
    }
  }
}

### ----------------------------------------------------------------------------
### .init.valid.Normult.Model
### @description    Initial validity check for parameters of a multivariate
###                 Normal mixture.
### @par    obj     a model object
### @return         An error if parameters fail certain conditions
### @detail         This initial validity check is called in the S4 constructor
###                 'model()' and ensures that the user constructs an inherently
###                 consistent model object.
###                 The parameter list must contain the foillowing elements:
###                     mu:     an r x K matrix containing 'numeric' or
###                             'integer' values
###                     sigma:  am r x r x K array containing 'numeric' or
###                             'integer' matrices, all symmetric/positive
###                             definite
### @see        ?model, ?vignette('finmix')
### @author     Lars Simon Zehnder
### ----------------------------------------------------------------------------
#' @noRd
".init.valid.Normult.Model" <- function(obj) {
  if (length(obj@par) > 0) {
    if (!"mu" %in% names(obj@par)) {
      stop(paste("Wrong specification of slot @par: ",
        "multivariate Normal mixtures need ",
        "a mean matrix named 'mu'.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!is.matrix(obj@par$mu)) {
      stop(paste("Wrong specification of slot @par: ",
        "mu is not a matrix. ",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!all(is.numeric(obj@par$mu) || is.numeric(obj@par$mu))) {
      stop(paste("Wrong specification of slot @par: ",
        "parameters must be of type 'numeric ",
        "or 'integer'.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!identical(dim(obj@par$mu), c(obj@r, obj@K))) {
      stop(paste("Wrong specification of slot @par: ",
        "mu must be a matrix of dimension r x K.",
        sep = ""
      ),
      call. = FALSE
      )
    }
    if (!"sigma" %in% names(obj@par)) {
      stop(paste("Wrong specification of slot @par: ",
        "multivariate Normal mixtures need ",
        "a variance-covariance array named ",
        "'sigma'",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!(is.numeric(obj@par$sigma) || is.integer(obj@par$mu))) {
      stop(paste("Wrong specification of slot @par: ",
        "parameters must be of type 'numeric' ",
        "or 'integer'.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!is.array(obj@par$sigma)) {
      stop(paste("Wrong specification of slot @par: ",
        "sigma is not an array.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!all(apply(obj@par$sigma, 3, isSymmetric))) {
      stop(paste("Wrong specification of slot @par: ",
        "sigma must contain K symmetric ",
        "r x r matrices.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!all(apply(obj@par$sigma, 3, function(x) {
      all(eigen(x)$values > 0)
    }))) {
      stop(paste("Wrong specification of slot @par: ",
        "sigma must contain K positive definite ",
        "r x r matrices.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!identical(dim(obj@par$sigma), c(obj@r, obj@r, obj@K))) {
      stop(paste("Wrong specification of slot @par: ",
        "sigma must be an array of dimension ",
        "r x r x K.",
        sep = ""
      ),
      call. = FALSE
      )
    }
  }
}

### ----------------------------------------------------------------------------
### .valid.Normult.Model
### @description    Initial validity check for parameters of a multivariate
###                 Normal mixture.
### @par    obj     a model object
### @return         An error if parameters fail necessary conditions and
###                 a warning if parameters fail consistency
### @detail         This validity check is called in the setters to ensure that
###                 slots can be changed without errors but help the user to
###                 end up with an inherently consistent model object.
###                 The parameter list must contain the foillowing elements:
###                     mu:     an r x K matrix containing 'numeric' or
###                             'integer' values
###                     sigma:  am r x r x K array containing 'numeric' or
###                             'integer' matrices, all symmetric/positive
###                             definite
### @see        ?model, ?vignette('finmix')
### @author     Lars Simon Zehnder
### ----------------------------------------------------------------------------
#' @noRd
".valid.Normult.Model" <- function(obj) {
  if (length(obj@par) > 0) {
    if (!"mu" %in% names(obj@par)) {
      warning(paste("Wrong specification of slot @par: ",
        "multivariate Normal mixtures need ",
        "a mean matrix named 'mu'.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!is.matrix(obj@par$mu)) {
      warning(paste("Wrong specification of slot @par: ",
        "mu is not a matrix. ",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!all(is.numeric(obj@par$mu) || is.numeric(obj@par$mu))) {
      stop(paste("Wrong specification of slot @par: ",
        "parameters must be of type 'numeric ",
        "or 'integer'.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!identical(dim(obj@par$mu), c(obj@r, obj@K))) {
      warning(paste("Wrong specification of slot @par: ",
        "mu must be a matrix of dimension r x K.",
        sep = ""
      ),
      call. = FALSE
      )
    }
    if (!"sigma" %in% names(obj@par)) {
      warning(paste("Wrong specification of slot @par: ",
        "multivariate Normal mixtures need ",
        "a variance-covariance array named ",
        "'sigma'",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!(is.numeric(obj@par$sigma) || is.integer(obj@par$mu))) {
      stop(paste("Wrong specification of slot @par: ",
        "parameters must be of type 'numeric' ",
        "or 'integer'.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!is.array(obj@par$sigma)) {
      warning(paste("Wrong specification of slot @par: ",
        "sigma is not an array.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!all(apply(obj@par$sigma, 3, isSymmetric))) {
      warning(paste("Wrong specification of slot @par: ",
        "sigma must contain K symmetric ",
        "r x r matrices.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!all(apply(obj@par$sigma, 3, function(x) {
      all(eigen(x)$values > 0)
    }))) {
      warning(paste("Wrong specification of slot @par: ",
        "sigma must contain K positive definite ",
        "r x r matrices.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!identical(dim(obj@par$sigma), c(obj@r, obj@r, obj@K))) {
      warning(paste("Wrong specification of slot @par: ",
        "sigma must be an array of dimension ",
        "r x r x K.",
        sep = ""
      ),
      call. = FALSE
      )
    }
  }
}

### ------------------------------------------------------------------------------
### .init.valid.Student.Model
### @description    Initial validity check for parameters of a univariate
###                 Student-t mixture.
### @par    obj     a model object
### @return         An error if parameters fail certain conditions
### @detail         This initial validity check is called in the S4 constructor
###                 'model()' and ensures that the user constructs an inherently
###                 consistent model object.
###                 The parameter list must contain the following elements:
###                     mu:     an 1 x K array, vector or matrix containing
###                             'numeric' or 'integer' values
###                     sigma:  an 1 x K array, vector or matrix containing
###                             'numeric' or 'integer' values, all positive
###                     df:     an 1 x K array, vector or matrix containing
###                             'numeric' or 'integer' values, all positive
### @see            ?model, ?vignette('finmix')
### @author         Lars Simon Zehnder
### -------------------------------------------------------------------------------
#' @noRd
".init.valid.Student.Model" <- function(obj) {
  if (length(obj@par) > 0) {
    if (!"mu" %in% names(obj@par)) {
      stop(paste("Wrong specification of slot @par: ",
        "univariate Normal mixtures need ",
        "a mean matrix named 'mu'.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!is.array(obj@par$mu) && !is.vector(obj@par$mu) &&
      !is.matrix(obj@par$mu)) {
      stop(paste("Wrong specification of slot @par: ",
        "mu must be either an array, a vector ",
        "or a matrix of dimension 1 x K. ",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!all(is.numeric(as.vector(obj@par$mu)) ||
      is.integer(as.vector(obj@par$mu)))) {
      stop(paste("Wrong specification of slot @par: ",
        "Parameters must be of type 'numeric' or ",
        "'integer'.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (length(obj@par$mu) != obj@K) {
      stop(paste("Wrong specification of slot @par: ",
        "mu must be a matrix of dimension 1 x K ",
        "or a vector of size K.",
        sep = ""
      ),
      call. = FALSE
      )
    }
    if (!"sigma" %in% names(obj@par)) {
      warning(paste("Wrong specification of slot @par: ",
        "univariate Student-t mixtures need ",
        "a variance vector named ",
        "'sigma'",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!all(is.numeric(as.vector(obj@par$sigma)) ||
      is.integer(as.vector(obj@par$sigma)))) {
      stop(paste("Wrong specification of slot @par: ",
        "Parameters must be of type 'numeric' or ",
        "'integer'.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (any(obj@par$sigma <= 0)) {
      stop(paste("Wrong specification of slot @par: ",
        "sigma must contain variances, all ",
        "positive.",
        sep = ""
      ),
      .call = FALSE
      )
    } else if (!is.array(obj@par$sigma) && is.vector(obj@par$sigma) &&
      is.matrix(obj@par$sigma)) {
      stop(paste("Wrong specification of slot @par: ",
        "sigma must be either an array, a vector, ",
        "or a matrix of dimension 1 x K.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (length(obj@par$sigma) != obj@K) {
      stop(paste("Wrong specification of slot @par: ",
        "sigma must be either an array, a vector, ",
        "or a matrix, ",
        "or a matrix of dimension ",
        "1 x K.",
        sep = ""
      ),
      call. = FALSE
      )
    }
    if (!"df" %in% names(obj@par)) {
      stop(paste("Wrong specification of slot @par: ",
        "Student-t mixtures need a degree of ",
        "freedom vector.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!all(is.numeric(as.vector(obj@par$df)) ||
      is.integer(as.vector(obj@par$df)))) {
      stop(paste("Wrong specification of slot @par: ",
        "Parameters must be of type 'numeric' or ",
        "'integer'.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (any(obj@par$df <= 0)) {
      stop(paste("Wrong specification of slot @par: ",
        "Degrees of freedom must be all positive.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (length(obj@par$df) != obj@K) {
      stop(paste("Wrong specification of slot @par: ",
        "df must be a vector or matrix of ",
        "dimension 1 x K",
        sep = ""
      ),
      call. = FALSE
      )
    }
  }
}

### ------------------------------------------------------------------------------
### .valid.Student.Model
### @description    Validity check for parameters of a univariate Student-t
###                 mixture.
### @par    obj     a model object
### @return         An error if parameters fail certain necessary conditions and
###                 a warning if parameters fail consistency.
### @detail         This validity check is called in the setters to ensure that
###                 slots can be changed without errors but help the user to
###                 end up with an inherently consistent model object.
###                 The parameter list must contain the following elements:
###                     mu:     an 1 x K array, vector or matrix containing
###                             'numeric' or 'integer' values
###                     sigma:  an 1 x K array, vector or matrix containing
###                             'numeric' or 'integer' values, all positive
###                     df:     an 1 x K array, vector or matrix containing
###                             'numeric' or 'integer' values, all positive
### @see            ?model, ?vignette('finmix')
### @author         Lars Simon Zehnder
### -------------------------------------------------------------------------------
#' @noRd
".valid.Student.Model" <- function(obj) {
  if (length(obj@par) > 0) {
    if (!"mu" %in% names(obj@par)) {
      warning(paste("Wrong specification of slot @par: ",
        "univariate Normal mixtures need ",
        "a mean matrix named 'mu'.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!is.array(obj@par$mu) && !is.vector(obj@par$mu) &&
      !is.matrix(obj@par$mu)) {
      warning(paste("Wrong specification of slot @par: ",
        "mu must be either an array, a vector ",
        "or a matrix of dimension 1 x K. ",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!all(is.numeric(as.vector(obj@par$mu)) ||
      is.integer(as.vector(obj@par$mu)))) {
      stop(paste("Wrong specification of slot @par: ",
        "Parameters must be of type 'numeric' or ",
        "'integer'.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (length(obj@par$mu) != obj@K) {
      warning(paste("Wrong specification of slot @par: ",
        "mu must be a matrix of dimension 1 x K ",
        "or a vector of size K.",
        sep = ""
      ),
      call. = FALSE
      )
    }
    if (!"sigma" %in% names(obj@par)) {
      warning(paste("Wrong specification of slot @par: ",
        "univariate Normal mictures need ",
        "a variance vector named ",
        "'sigma'",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!all(is.numeric(as.vector(obj@par$sigma)) ||
      is.integer(as.vector(obj@par$sigma)))) {
      stop(paste("Wrong specification of slot @par: ",
        "Parameters must be of type 'numeric' or ",
        "'integer'.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (any(obj@par$sigma <= 0)) {
      stop(paste("Wrong specification of slot @par: ",
        "sigma must contain variances, all ",
        "positive.",
        sep = ""
      ),
      .call = FALSE
      )
    } else if (is.array(obj@par$sigma) && is.vector(obj@par$sigma) &&
      is.matrix(obj@par$sigma)) {
      warning(paste("Wrong specification of slot @par: ",
        "sigma must be either an array, a vector, ",
        "or a matrix of dimension 1 x K.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (length(obj@par$sigma) != obj@K) {
      warning(paste("Wrong specification of slot @par: ",
        "sigma must be either an array, a vector, ",
        "or a matrix, ",
        "or a matrix of dimension ",
        "1 x K.",
        sep = ""
      ),
      call. = FALSE
      )
    }
    if (!"df" %in% names(obj@par)) {
      warning(paste("Wrong specification of slot @par: ",
        "Student-t mixtures need a degree of ",
        "freedom vector.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!all(is.numeric(as.vector(obj@par$df)) ||
      is.integer(as.vector(obj@par$df)))) {
      stop(paste("Wrong specification of slot @par: ",
        "Parameters must be of type 'numeric' or ",
        "'integer'.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (any(obj@par$df <= 0)) {
      stop(paste("Wrong specification of slot @par: ",
        "Degrees of freedom must be all positive.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (length(obj@par$df) != obj@K) {
      warning(paste("Wrong specification of slot @par: ",
        "df must be a vector or matrix of ",
        "dimension 1 x K",
        sep = ""
      ),
      call. = FALSE
      )
    }
  }
}

### ----------------------------------------------------------------------------
### .init.valid.Studmult.Model
### @description    Initial validity check for parameters of a multivariate
###                 Student-t mixture.
### @par    obj     a model object
### @return         An error if parameters fail certain conditions
### @detail         This initial validity check is called in the S4 constructor
###                 'model()' and ensures that the user constructs an inherently
###                 consistent model object.
###                 The parameter list must contain the foillowing elements:
###                     mu:     an r x K matrix containing 'numeric' or
###                             'integer' values
###                     sigma:  an r x r x K array containing 'numeric' or
###                             'integer' matrices, all symmetric/positive
###                             definite
###                     df:     an 1 x K array, vector or matrix containing
###                             'numeric' or 'integer', all positive
### @see        ?model, ?vignette('finmix')
### @author     Lars Simon Zehnder
### ----------------------------------------------------------------------------
#' @noRd
".init.valid.Studmult.Model" <- function(obj) {
  if (length(obj@par) > 0) {
    if (!"mu" %in% names(obj@par)) {
      stop(paste("Wrong specification of slot @par: ",
        "multivariate Student-t mixtures need ",
        "a mean matrix named 'mu'.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!is.matrix(obj@par$mu)) {
      stop(paste("Wrong specification of slot @par: ",
        "mu is not a matrix. ",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!all(is.numeric(obj@par$mu) || is.numeric(obj@par$mu))) {
      stop(paste("Wrong specification of slot @par: ",
        "parameters must be of type 'numeric ",
        "or 'integer'.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!identical(dim(obj@par$mu), c(obj@r, obj@K))) {
      stop(paste("Wrong specification of slot @par: ",
        "mu must be a matrix of dimension r x K.",
        sep = ""
      ),
      call. = FALSE
      )
    }
    if (!"sigma" %in% names(obj@par)) {
      stop(paste("Wrong specification of slot @par: ",
        "multivariate Student-t mictures need ",
        "a variance-covariance array named ",
        "'sigma'",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!(is.numeric(obj@par$sigma) || is.integer(obj@par$mu))) {
      stop(paste("Wrong specification of slot @par: ",
        "parameters must be of type 'numeric' ",
        "or 'integer'.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!is.array(obj@par$sigma)) {
      stop(paste("Wrong specification of slot @par: ",
        "sigma is not an array.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!all(apply(obj@par$sigma, 3, isSymmetric))) {
      stop(paste("Wrong specification of slot @par: ",
        "sigma must contain K symmetric ",
        "r x r matrices.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!all(apply(obj@par$sigma, 3, function(x) {
      all(eigen(x)$values > 0)
    }))) {
      stop(paste("Wrong specification of slot @par: ",
        "sigma must contain K positive definite ",
        "r x r matrices.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!identical(dim(obj@par$sigma), c(obj@r, obj@r, obj@K))) {
      stop(paste("Wrong specification of slot @par: ",
        "sigma must be an array of dimension ",
        "r x r x K.",
        sep = ""
      ),
      call. = FALSE
      )
    }
    if (!"df" %in% names(obj@par)) {
      warning(paste("Wrong specification of slot @par: ",
        "Student-t mixtures need a degree of ",
        "freedom vector.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!all(is.numeric(as.vector(obj@par$df)) ||
      is.integer(as.vector(obj@par$df)))) {
      stop(paste("Wrong specification of slot @par: ",
        "Parameters must be of type 'numeric' or ",
        "'integer'.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (any(obj@par$df <= 0)) {
      stop(paste("Wrong specification of slot @par: ",
        "Degrees of freedom must be all positive.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (length(obj@par$df) != obj@K) {
      warning(paste("Wrong specification of slot @par: ",
        "df must be a vector or matrix of ",
        "dimension 1 x K",
        sep = ""
      ),
      call. = FALSE
      )
    }
  }
}

### ----------------------------------------------------------------------------
### .valid.Studmult.Model
### @description    Initial validity check for parameters of a multivariate
###                 Student-t mixture.
### @par    obj     a model object
### @return         An error if parameters fail necessary conditions and
###                 a warning if parameters fail consistency
### @detail         This validity check is called in the setters to ensure that
###                 slots can be changed without errors but help the user to
###                 end up with an inherently consistent model object.
###                 The parameter list must contain the foillowing elements:
###                     mu:     an r x K matrix containing 'numeric' or
###                             'integer' values
###                     sigma:  am r x r x K array containing 'numeric' or
###                             'integer' matrices, all symmetric/positive
###                             definite
###                     df:     an 1 x K array, vector or matrix containing
###                             'numeric' or 'integer' values, all positive
### @see        ?model, ?vignette('finmix')
### @author     Lars Simon Zehnder
### ----------------------------------------------------------------------------
#' @noRd
".valid.Studmult.Model" <- function(obj) {
  if (length(obj@par) > 0) {
    if (!"mu" %in% names(obj@par)) {
      warning(paste("Wrong specification of slot @par: ",
        "multivariate Student-t mixtures need ",
        "a mean matrix named 'mu'.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!is.matrix(obj@par$mu)) {
      warning(paste("Wrong specification of slot @par: ",
        "mu is not a matrix. ",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!all(is.numeric(obj@par$mu) || is.numeric(obj@par$mu))) {
      stop(paste("Wrong specification of slot @par: ",
        "parameters must be of type 'numeric ",
        "or 'integer'.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!identical(dim(obj@par$mu), c(obj@r, obj@K))) {
      warning(paste("Wrong specification of slot @par: ",
        "mu must be a matrix of dimension r x K.",
        sep = ""
      ),
      call. = FALSE
      )
    }
    if (!"sigma" %in% names(obj@par)) {
      warning(paste("Wrong specification of slot @par: ",
        "multivariate Student-t mictures need ",
        "a variance-covariance array named ",
        "'sigma'",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!(is.numeric(obj@par$sigma) || is.integer(obj@par$mu))) {
      stop(paste("Wrong specification of slot @par: ",
        "parameters must be of type 'numeric' ",
        "or 'integer'.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!is.array(obj@par$sigma)) {
      warning(paste("Wrong specification of slot @par: ",
        "sigma is not an array.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!all(apply(obj@par$sigma, 3, isSymmetric))) {
      warning(paste("Wrong specification of slot @par: ",
        "sigma must contain K symmetric ",
        "r x r matrices.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!all(apply(obj@par$sigma, 3, function(x) {
      all(eigen(x)$values > 0)
    }))) {
      warning(paste("Wrong specification of slot @par: ",
        "sigma must contain K positive definite ",
        "r x r matrices.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!identical(dim(obj@par$sigma), c(obj@r, obj@r, obj@K))) {
      warning(paste("Wrong specification of slot @par: ",
        "sigma must be an array of dimension ",
        "r x r x K.",
        sep = ""
      ),
      call. = FALSE
      )
    }
    if (!"df" %in% names(obj@par)) {
      warning(paste("Wrong specification of slot @par: ",
        "Student-t mixtures need a degree of ",
        "freedom vector.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (!all(is.numeric(as.vector(obj@par$df)) ||
      is.integer(as.vector(obj@par$df)))) {
      stop(paste("Wrong specification of slot @par: ",
        "Parameters must be of type 'numeric' or ",
        "'integer'.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (any(obj@par$df <= 0)) {
      stop(paste("Wrong specification of slot @par: ",
        "Degrees of freedom must be all positive.",
        sep = ""
      ),
      call. = FALSE
      )
    } else if (length(obj@par$df) != obj@K) {
      warning(paste("Wrong specification of slot @par: ",
        "df must be a vector or matrix of ",
        "dimension 1 x K",
        sep = ""
      ),
      call. = FALSE
      )
    }
  }
}

### Additional functions
#' Returns all univariate distributions
#' 
#' @description 
#' For internal usage only. 
#' 
#' @return A character vector containing all univariate distributions.
#' @noRd
".get.univ.Model" <- function() {
  univ <- c(
    "poisson", "cond.poisson",
    "binomial", "exponential",
    "normal", "student"
  )
  return(univ)
}

#' Returns all multivariate distributions
#' 
#' @description 
#' For internal usage only. 
#' 
#' @return A character vector containing all multivariate distributions.
#' @noRd
".get.multiv.Model" <- function() {
  multiv <- c("normult", "studmult")
  return(multiv)
}
