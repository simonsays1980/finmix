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

# CLASSIFICATION OF OBSEVRATIONS USING FULLY SPECIFIED MODEL
#' Finmix `dataclass` class
#' 
#' @description
#' Stores objects to classify observations using a fully specified mixture 
#' model. If the indicators a finite mixture model is fully specified as then 
#' the likelihood can be calculated for each observation depending on the 
#' component it stems from. 
#' 
#' @slot logpy An array containing the logarithmized 
#' @slot prob An array storing the probability classification matrix that 
#'   defines for each observation the probability of belonging to component 
#'   `k`. Henceforth, each row sums to one. The dimension of this array is 
#'   `N x K`. 
#' @slot mixlik A numeric storing the logarithm of the mixture likelihood 
#'   evaluated at certain parameters `par` from a finmix `model` object and 
#'   corresponding `weights`.
#' @slot entropy A numeric storing the entropy of the classification.
#' @slot loglikcd An array storing the logarithm of the conditional likelihood 
#'   of each component parameter, if indicators have not been simulated. The 
#'   array has dimension `1 x K`.
#' @slot postS A numeric storing the posterior probability of the indicators 
#'   `S` in the data, if indicators have been simulated.
#' @exportClass dataclass
#' @name dataclass_class
#' 
#' @seealso 
#' * [fdata][fdata_class] for the class holding the data
#' * [model][model_class] for the class defining a finite mixture model
#' * [dataclass()] for the constructor of this class
#' 
#' @references 
#' Frühwirth-Schnatter, S. (2006), "Finite Mixture and Markov Switching Models"
.dataclass <- setClass("dataclass",
  representation(
    logpy     = "array",
    prob      = "array",
    mixlik    = "numeric",
    entropy   = "numeric",
    loglikcd  = "array",
    postS     = "numeric"
  ),
  validity = function(object) {
    ## else: OK
    TRUE
  },
  prototype(
    logpy = array(),
    prob = array(),
    mixlik = numeric(),
    entropy = numeric(),
    loglikcd = array(),
    postS = numeric()
  )
)

#' Finmix `dataclass` constructor
#' 
#' Calling [dataclass()] classifies data using a fully specified mixture model. 
#' Henceforth, the finite mixture model `model` must be fully specified, i.e. 
#' containing parameters in slot `@@par`, weights in slot `@@weight` and 
#' indicators in slot `@@S` of the corresponding `fdata` object.
#' 
#' @param fdata An `fdata` object containing observations in slot `@@y` and 
#'   indicators in slot `@@S`.
#' @param model A `model` object containing parameters in slot  `@@par` and 
#'   and weights in slot `@@weight`.
#' @param simS A logical defining, if the indicators `S` should be simulated.
#' @return A `dataclass` object containing the classification matrix, 
#'   model log-likelihood, entropy and indicators, if the latter have been 
#'   simulated.
#' @export
#' 
#' @seealso 
#' * [dataclass][dataclass_class] for the class definition
#' 
#' #' @references 
#' Frühwirth-Schnatter, S. (2006), "Finite Mixture and Markov Switching Models"
"dataclass" <- function(fdata = NULL, model = NULL, simS = FALSE) {
  .check.fdata.model.Dataclass(fdata, model)
  .check.model.Dataclass(model)
  .valid.fdata.model.Prior(fdata, model)
  if (hasS(fdata)) {
    classm <- getColS(fdata)
  }
  datam <- getColY(fdata)
  K <- model@K
  dist <- model@dist
  lik.list <- .liklist.Dataclass(fdata, model)
  ## The following is only used as a temporary solution as long
  ## as the 'Markov' model for the indicators is not yet implemented.
  ## Check attribute 'indicmod' in 'model' argument
  if (model@indicmod != "multinomial") {
    model@indicmod <- "multinomial"
  }
  if (!model@indicfix) {
    if (model@indicmod == "multinomial") {
      .multinomial.Dataclass(
        fdata, model, lik.list,
        simS
      )
    } ## else: implemented later: Markov model for S
  } else { ## indicfix == TRUE
    .indicfix.Dataclass(fdata, model, lik.list)
  }
}

#' Shows a summary of a `dataclass` object.
#' 
#' Calling [show()] on a `dataclass` object gives an overview 
#' of the slots of this class.
#' 
#' @param object A `dataclass` object.
#' @returns A console output listing the slots and summary information about
#'   each of them. 
#' @exportMethod  show
#' @describeIn dataclass_class
setMethod(
  "show", "dataclass",
  function(object) {
    has.loglikcd <- !all(is.na(object@loglikcd))
    has.mixlik <- !all(is.na(object@mixlik))
    has.entropy <- !all(is.na(object@entropy))
    has.postS <- !all(is.na(object@postS))
    cat("Object 'dataclass'\n")
    cat(
      "     logpy       :",
      paste(dim(object@logpy), collapse = "x"),
      "\n"
    )
    cat(
      "     prob        :",
      paste(dim(object@prob), collapse = "x"),
      "\n"
    )
    if (has.mixlik) {
      cat("     mixlik      :", object@mixlik, "\n")
    }
    if (has.entropy) {
      cat("     entropy     :", object@entropy, "\n")
    }
    if (has.loglikcd) {
      cat(
        "     loglikcd    :",
        paste(dim(object@loglikcd), collapse = "x"),
        "\n"
      )
    }
    if (has.postS) {
      cat("     postS       :", object@postS, "\n")
    }
  }
)

## Getters ##
#' Getter method of `dataclass` class.
#' 
#' Returns the `logpy` slot.
#' 
#' @param object An `dataclass` object.
#' @returns The `logpy` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' # Generate a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 0.7)), K = 2)
#' # Simulate data from the model.
#' f_data <- simulate(f_model)
#' # Classify observations
#' f_dataclass <- dataclass(f_data, f_model, simS = FALSE)
#' getLogpy(f_datamoms)
#' 
#' @seealso 
#' * [dataclass][dataclass_class] for the base class
#' * [dataclass()] for the constructor of the `dataclass` class
setMethod(
  "getLogpy", "dataclass",
  function(object) {
    return(object@logpy)
  }
)

#' Getter method of `dataclass` class.
#' 
#' Returns the `prob` slot.
#' 
#' @param object An `dataclass` object.
#' @returns The `prob` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' # Generate a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 0.7)), K = 2)
#' # Simulate data from the model.
#' f_data <- simulate(f_model)
#' # Classify observations
#' f_dataclass <- dataclass(f_data, f_model, simS = FALSE)
#' getProb(f_datamoms)
#' 
#' @seealso 
#' * [dataclass][dataclass_class] for the base class
#' * [dataclass()] for the constructor of the `dataclass` class
setMethod(
  "getProb", "dataclass",
  function(object) {
    return(object@prob)
  }
)

#' Getter method of `dataclass` class.
#' 
#' Returns the `mixlik` slot.
#' 
#' @param object An `dataclass` object.
#' @returns The `mixlik` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' # Generate a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 0.7)), K = 2)
#' # Simulate data from the model.
#' f_data <- simulate(f_model)
#' # Classify observations
#' f_dataclass <- dataclass(f_data, f_model, simS = FALSE)
#' getMixlik(f_datamoms)
#' 
#' @seealso 
#' * [dataclass][dataclass_class] for the base class
#' * [dataclass()] for the constructor of the `dataclass` class<
setMethod(
  "getMixlik", "dataclass",
  function(object) {
    return(object@mixlik)
  }
)

#' Getter method of `dataclass` class.
#' 
#' Returns the `entropy` slot.
#' 
#' @param object An `dataclass` object.
#' @returns The `entropy` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' # Generate a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 0.7)), K = 2)
#' # Simulate data from the model.
#' f_data <- simulate(f_model)
#' # Classify observations
#' f_dataclass <- dataclass(f_data, f_model, simS = FALSE)
#' getEntropy(f_datamoms)
#' 
#' @seealso 
#' * [dataclass][dataclass_class] for the base class
#' * [dataclass()] for the constructor of the `dataclass` class
setMethod(
  "getEntropy", "dataclass",
  function(object) {
    return(object@entropy)
  }
)

#' Getter method of `dataclass` class.
#' 
#' Returns the `loglikcd` slot. Note that this slot is only non-null, if the 
#' indicators have not been simulated.
#' 
#' @param object An `dataclass` object.
#' @returns The `loglikcd` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' # Generate a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 0.7)), K = 2)
#' # Simulate data from the model.
#' f_data <- simulate(f_model)
#' # Classify observations
#' f_dataclass <- dataclass(f_data, f_model, simS = FALSE)
#' getLoglikcd(f_datamoms)
#' 
#' @seealso 
#' * [dataclass][dataclass_class] for the base class
#' * [dataclass()] for the constructor of the `dataclass` class
setMethod(
  "getLoglikcd", "dataclass",
  function(object) {
    return(object@loglikcd)
  }
)

#' Getter method of `dataclass` class.
#' 
#' Returns the `postS` slot. Note that this slot is only non-null, if the 
#' indicators have been simulated.
#' 
#' @param object An `dataclass` object.
#' @returns The `postS` slot of the `object`.
#' @noRd
#' 
#' @examples 
#' # Generate a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 0.7)), K = 2)
#' # Simulate data from the model.
#' f_data <- simulate(f_model)
#' # Classify observations
#' f_dataclass <- dataclass(f_data, f_model, simS = TRUE)[[1]]
#' getPostS(f_datamoms)
#' 
#' @seealso 
#' * [dataclass][dataclass_class] for the base class
#' * [dataclass()] for the constructor of the `dataclass` class
setMethod(
  "getPostS", "dataclass",
  function(object) {
    return(object@postS)
  }
)

## No setters as users are not intended to manipulate ##
## this object.                                       ##

### Private functions
### These functions are not exported.

### Checking
### Check fdata/model: 'fdata' must be an object of class
### 'fdata'. Further, this object must be valid and must
### contain data in @y. The 'fdata' object and the 'model'
### object must be consistent to each other, i.e. the 'model'
### object must have defined a distribution in @dist that
### conforms with the dimension @r of the #fdata' object.
#' Checking `fdata` object and `model` object for `dataclass`
#' 
#' For internal usage only. This function checks an `fdata` object and a 
#' `model` object in regard to consistency. First of all the data dimensions 
#' must fit between the two object, meaning that if `@@r>1` in the `fdata` 
#' object the model object must possess a `@@dist` slot with an appropriate 
#' distribution for multivariate data. The `fdata` object must contain data in 
#' its slot `@@y`. As a first safeguard this function checks if the first 
#' argument is indeed an `fdata` object. 
#' 
#' @param fdata.obj An `fdata` object. 
#' @param model.obj A `model` object containing a specified distribution, 
#'   parameters and weigths.
#' @return None. If the checks do not run through, an error is thrown.
#' @noRd
#' 
#' @seealso 
#' * [fdata][fdata_class] for the `fdata` class definition
#' * [model][model_class] for the `model` class definition
#' * [dataclass()] for the calling function
".check.fdata.model.Dataclass" <- function(fdata.obj, model.obj) {
  if (class(fdata.obj) != "fdata") {
    stop(paste("Wrong argument in 'dataclass()'. First ",
      "argument must be an object of class 'fdata'",
      sep = ""
    ))
  }
  .valid.Fdata(fdata.obj)
  hasY(fdata.obj, verbose = TRUE)
}

### Check model: 'model' must be an object of class 'model'.
### Furthermore, it must be valid and contain specified
### parameters in @par and weights in @weight.
#' Checking `model` object for `dataclass` constructor 
#' 
#' For internal usage only. This function checks if the `model` object passed
#' in by the user is first of all indeed a finmix `model` object. Furthermore, 
#' it is checked if the model is fully specified meaning that parameters are 
#' defined in slot `@@par` and weights in slot `@@weight`. 
#' 
#' @param model.obj A `model` object. Must be fully specified.
#' @return None. If the checks do not pass, an error is thrown. 
#' @noRd
#' 
#' @seealso 
#' * [fdata][fdata_class] for the `fdata` class definition
#' * [model][model_class] for the `model` class definition
#' * [dataclass()] for the calling function
".check.model.Dataclass" <- function(model.obj) {
  if (class(model.obj) != "model") {
    stop(paste("Wrong argument in 'dataclass()'. Second ",
      "argument must be an object of class 'model'.",
      sep = ""
    ))
  }
  .valid.Model(model.obj)
  hasPar(model.obj, verbose = TRUE)
  hasWeight(model.obj, verbose = TRUE)
}

### Check indicators: Indicators must have as many different factors
### as @K in the 'model' object. Further, values must be out of the
### sequence 1, ..., K.
#' Checking indicators for `dataclass` constructor 
#' 
#' For internal usage only. This function checks if the indicators stored in 
#' the `fdata` object are correctly specified meaning if the indicator values 
#' are indeed from as many components as specifed in the slot `@@K` of the 
#' corresponding model object.  
#' 
#' @param fdata.obj An `fdata` object containing the indicators in its slot 
#'   `@@S`.
#' @param model.obj A `model` object. Must be fully specified. 
#' @return None. If the checks do not pass, an error is thrown. 
#' @noRd
#' 
#' @seealso 
#' * [fdata][fdata_class] for the `fdata` class definition
#' * [model][model_class] for the `model` class definition
#' * [dataclass()] for the calling function
".check.S.Dataclass" <- function(fdata.obj, model.obj) {
  values <- levels(as.factor(fdata.obj@S))
  if (!identical(range(fdata.obj@S), range(seq(1, model.obj@K)))) {
    stop(paste("Wrong specification of slot 'S' in 'fdata' ",
      "object. Indicators must have integer values ",
      "in the range 1 to slot 'K' of 'model' object.",
      sep = ""
    ))
  }
}

#' Checking Student-t and normal `model` objects for `dataclass` constructor 
#' 
#' For internal usage only. Thiss function checks if the `model` object passed
#' in by the user is correctly specified in case the slot `@@dist` is one of 
#' `normult` or `studmult`. Correctly specified for data classification means 
#' that the slots `@@sigmainv` and `@@logdet` are non-null. `@@sigmainv` is the 
#' inverse of the variance-covariance matrix of a multivariate normal or 
#' Student-t distribution. Slot `@@logdet` defines the logarithm of the 
#' determinant of the inverse of the variance-covariance matrix. If these slots 
#' are not specified this function specifies these slots for the user.
#' 
#' @param model.obj A `model` object. Must be fully specified.
#' @return The passed-in `model` object by the user possibly enriched by slots 
#'   `@@sigmainv` and `@@logdet`. If the checks do not pass, an error is thrown. 
#' @noRd
#' 
#' @seealso 
#' * [fdata][fdata_class] for the `fdata` class definition
#' * [model][model_class] for the `model` class definition
#' * [dataclass()] for the calling function
".check.Logdet.Norstud" <- function(model.obj) {
  has.sigmainv <- "sigmainv" %in% names(model.obj@par)
  has.logdet   <- "logdet" %in% names(model.obj@par)
  r            <- model.obj@r 
  K            <- model.obj@K
  if (has.sigmainv && has.logdet) {
    return(model.obj)
  } else {
    qinv <- array(0, dim = c(r, r, K))
    logdetq <- array(0, dim = c(1, K))
    for (k in 1:K) {
      qinv[, , k] <- solve(model@par$sigma[, , k])
      logdetq[k] <- log(det(qinv[, , k]))
    }
    model.obj@par$qinv <- qinv
    model.obj@par$logdetq <- logdetq
    return(model.obj)
  }
}

### Logic
### Logic liklist
### Compute the likelihood l(y_i|theta_k) for all i and k
### lik.list is a 'list' object containing ##
###      'lh'          exp(llh - maxl), an N x K  'matrix'
###      'maxl'        the maximum likelihood, an N x 1 'vector'
###      'llh'         the likelihood, a N x K 'matrix'
#' Computes the log-likelihood for data classification
#' 
#' @description
#' For internal usage only. This function calls the appropriate function for 
#' each finite mixture model specified in `model.obj`. 
#' 
#' @param fdata.obj An `fdata` object with non-empty slot `@@y`.
#' @param model.obj A `model` object. Must be fully specified.
#' @return A list containing the likelihood, the maximum likelihood and the 
#'   log-likelihood.
#' @noRd
#' 
#' @seealso 
#' * [fdata][fdata_class] for the `fdata` class definition
#' * [model][model_class] for the `model` class definition
#' * [dataclass()] for the calling function
".liklist.Dataclass" <- function(fdata.obj, model.obj) {
  K <- model.obj@K
  N <- fdata.obj@N
  dist <- model.obj@dist
  datam <- getColY(fdata.obj)
  if (dist == "normal") {
    .likelihood.normal(
      datam, model.obj@par$mu,
      model.obj@par$sigma
    )
  } else if (dist == "student") {
    .likelihood.student(
      datam, model.obj@par$mu,
      model.obj@par$sigma,
      model.obj@par$df
    )
  } else if (dist == "exponential") {
    .likelihood.exponential(datam, model.obj@par$lambda)
  } else if (dist == "poisson" || dist == "cond.poisson") {
    ## should give a N x K 'matrix' object
    lambda <- model.obj@par$lambda
    if (hasExp(fdata.obj)) {
      expos <- getExp(fdata.obj)
      lambda <- matrix(lambda,
        nrow = N,
        ncol = K, byrow = TRUE
      )
      lambda <- apply(lambda, 2, "*", expos)
    } else {
      lambda <- matrix(lambda,
        nrow = 1,
        ncol = K
      )
    }
    .likelihood.poisson(datam, lambda)
  } else if (dist == "binomial") {
    .likelihood.binomial(datam, fdata.obj@T, model.obj@par$p)
  } else if (dist == "normult") {
    model.obj <- .check.Logdet.Norstud(model.obj)
    .likelihood.normult(
      datam, model.obj@par$mu,
      model.obj@par$sigmainv,
      model.obj@par$logdet
    )
  } else if (dist == "studmult") {
    model.obj <- .check.Logdet.Norstud(model.obj)
    .likelihood.studmult(
      datam, model.obj@par$mu,
      model.obj@par$sigmainv,
      model.obj@par$logdet,
      model.obj@par$df
    )
  }
}

#' Data classification for a multinomial indicator model
#' 
#' @description
#' For internal usage only. This function constructs the `dataclass` object 
#' and, if specified, simulates indicators `S`. A corresponding classification 
#' for a Markov indicator model is not (yet) implemented.
#' 
#' @param fdata.obj An `fdata` object with non-empty slot `@@y`.
#' @param model.obj A `model` object. Must be fully specified.
#' @param lik.list A list containing the likelihood, maximum likelihood and 
#'   log-likelihood for the data using the specified model. 
#' @param A logical specifying, if indicators should be simulated.
#' @return A `dataclass` object containing the classifications of the data, if 
#'   `simS` was set to `TRUE`, and the likelihood values.
#' @noRd
#' 
#' @seealso 
#' * [fdata][fdata_class] for the `fdata` class definition
#' * [model][model_class] for the `model` class definition
#' * [dataclass()] for the calling function
".multinomial.Dataclass" <- function(fdata.obj, model.obj,
                                     lik.list, simS) {
  N <- fdata.obj@N
  K <- model.obj@K
  mixlik.list <- .mixlik.Dataclass(model.obj, lik.list,
    prob = TRUE
  )
  mixlik <- mixlik.list$mixlik
  p <- mixlik.list$p
  if (simS && K > 1) {
    sim.S <- .simulate.S.Dataclass(p, K, N)
    S <- sim.S$S
    postS <- sim.S$postS
  } else {
    ## compute complete data likelihood in case
    ## indicators were not simulated
    if (hasS(fdata.obj) && K > 1) {
      classm <- getColS(fdata.obj)
      .check.S.Dataclass(fdata.obj, model.obj)
      loglikcd <- matrix(0, nrow = 1, ncol = K)
      for (k in seq(1, K)) {
        loglikcd[k] <- sum(lik.list$llh[classm == k, k])
      }
    } else { ## no indicators given or no mixture ##
      loglikcd <- matrix(mixlik)
    }
  }
  ## compute entropy
  logp <- matrix(0, nrow = N, ncol = K)
  for (k in seq(1, K)) {
    logp[p[, k] == 0, k] <- -99
    logp[p[, k] != 0, k] <- log(p[p[, k] != 0, k])
  }
  entropy <- (-1) * sum(logp * p)
  if (simS) {
    datac.obj <- .dataclass(
      logpy = lik.list$llh,
      prob = p, mixlik = mixlik,
      entropy = entropy,
      loglikcd = matrix(), postS = postS
    )
    l <- list(dataclass = datac.obj, S = as.integer(S))
    return(l)
  } else {
    .dataclass(
      logpy = lik.list$llh,
      prob = p, mixlik = mixlik,
      entropy = entropy,
      loglikcd = loglikcd, postS = numeric()
    )
  }
}

#' Computes the mixture likelihood for data classification
#' 
#' @description
#' For internal usage only. This function computes the mixture likelihood for 
#' the finite mixture model specified in `model.obj` using the likelihoods 
#' of each single component and the weights specified in slot `@@weight` of the 
#' `model` object.  
#' 
#' @param model.obj A `model` object. Must be fully specified.
#' @param lik.list A list containing the likelihood, maximum likelihood and the 
#'   log-likelihood of the data for each component of the finite mixture.
#' @param prob A logical indicating, if the probability classification matrix 
#'   should be computed.
#' @return A matrix of dimensions `N x 1` containing the mixture likelihood 
#'   for each observation. If `prob` ws set to `TRUE`, a list is returned 
#'   containing the mixture likelihood and in addition the probability 
#'   classification matrix of dimension `N x K`.
#' @noRd
#' 
#' @seealso 
#' * [fdata][fdata_class] for the `fdata` class definition
#' * [model][model_class] for the `model` class definition
#' * [dataclass()] for the calling function
".mixlik.Dataclass" <- function(model.obj, lik.list, prob = FALSE) {
  ## p is an N x K matrix ##
  p <- t(apply(lik.list$lh, 1, "*", model.obj@weight))
  ## sump is an N x 1 matrix ##
  sump <- apply(p, 1, sum)
  ## lsump is an N x 1 matrix ##
  lsump <- log(sump) + lik.list$maxl
  mixlik <- sum(lsump) ## numeric
  if (prob) {
    ## p is the N x K probability classification matrix ##
    p <- apply(p, 2, "/", sump)
    return(list(mixlik = mixlik, p = p))
  } else {
    return(mixlik)
  }
}

#' Simulate classification for a finite mixture model
#' 
#' @description
#' For internal usage only. This function simulates the indicators for a finite 
#' mixture model using the probability classification matrix. 
#' 
#' @param p A matrix containing the classification probabilities for each 
#'   component `K` of the finite mixture and each observation. Dimension is 
#'   `N x K`.
#' @param K A numeric specifying the value range for the simulated indicators. 
#'   Corresponds to the number of components in the finite mixture model.
#' @param N A numeric specifying the number of indicators to be simulated. 
#' @return A list containing the simulated indicators together with the 
#'   posterior log-likelihood of the simulated indicators.
#' @noRd
#' @importFrom stats runif
#' @seealso 
#' * [fdata][fdata_class] for the `fdata` class definition
#' * [model][model_class] for the `model` class definition
#' * [dataclass()] for the calling function
".simulate.S.Dataclass" <- function(p, K, N) {
  ## Simulate classifications from classification probability
  ## matrix
  rnd <- runif(N)
  S <- t(apply(p, 1, cumsum)) < matrix(rnd,
    nrow = length(rnd),
    ncol = K
  )
  S <- matrix(apply(S, 1, sum)) + 1
  Sm <- matrix(S, nrow = nrow(S), ncol = K)
  Compm <- matrix(seq(1, K), nrow = nrow(S), ncol = K, byrow = TRUE)
  postS <- sum(log(apply((Sm == Compm) * p, 1, sum)))
  sim.S <- list(S = S, postS = postS)
  return(sim.S)
}

#' Classification with fixed indicators
#' 
#' @description
#' For internal usage only. This function classifies data from a finite 
#' mixture model with fixed indicators.
#' 
#' @param fdata.obj An `fdata` object with non-empty slot `@@y`.
#' @param model.obj A `model` object. Must be fully specified. The slot 
#'   `@@indicfix` must be `TRUE`.
#' @param lik.list A list containing the likelihood, maximum likelihood, and 
#'   log-likelihood for the data in the `fdata` object.
#' @return An object of class `dataclass` containing the likelihood values for 
#'   the finite mixture model and the the data
#' @noRd
#' 
#' @seealso 
#' * [fdata][fdata_class] for the `fdata` class definition
#' * [model][model_class] for the `model` class definition
#' * [dataclass()] for the calling function
".indicfix.Dataclass" <- function(fdata.obj, model.obj, lik.list) {
  K <- model.obj@K
  mixlik <- .mixlik.Dataclass(model.obj, lik.list)
  if (hasS(fdata.obj) && K > 1) {
    .check.S.Dataclass(fdata.obj, model.obj)
    classm <- getColS(fdata.obj)
    loglikcd <- matrix(0, nrow = 1, ncol = K)
    for (k in seq(1, K)) {
      loglikcd[k] <- sum(lik.list$llh[classm == k, k])
    }
  } else { ## no indicators given or no mixture ##
    loglikcd <- matrix(mixlik)
  }
  .dataclass(
    logpy = lik.list$llh,
    prob = matrix(), mixlik = numeric(),
    entropy = numeric(),
    loglikcd = loglikcd, postS = numeric()
  )
}
