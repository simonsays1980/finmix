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
setMethod(
  "getLogpy", "dataclass",
  function(object) {
    return(object@logpy)
  }
)
setMethod(
  "getProb", "dataclass",
  function(object) {
    return(object@prob)
  }
)
setMethod(
  "getMixlik", "dataclass",
  function(object) {
    return(object@mixlik)
  }
)
setMethod(
  "getEntropy", "dataclass",
  function(object) {
    return(object@entropy)
  }
)
setMethod(
  "getLoglikcd", "dataclass",
  function(object) {
    return(object@loglikcd)
  }
)
setMethod(
  "getPostS", "dataclass",
  function(object) {
    return(object@postS)
  }
)

## No setters as users are not intended to mnaipulate ##
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

".check.Logdet.Norstud" <- function(model.obj) {
  has.sigmainv <- "sigmainv" %in% names(model.obj@par)
  has.logdet <- "logdet" %in% names(model.obj@par)
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
###      'maxl'        the maximum likelihood, an 1 x K 'vector'
###      'llh'         the likelihood, a N x K 'matrix'
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
