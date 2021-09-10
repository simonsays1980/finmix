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

"mixturemcmc" <- function(fdata, model, prior, mcmc) {
  ## Check arguments
  mcmc <- .check.args.Mixturemcmc(fdata, model, prior, mcmc, nargs())

  ## Default ordering for MCMC: bycolumn
  setBycolumn(fdata) <- TRUE
  ######################### MCMC SAMPLING #############################
  ## Set the indicators as a default to one for K == 1
  if (model@K == 1) {
    fdata@S <- matrix(1, nrow = fdata@N, ncol = 1)
  }
  dist <- model@dist
  if (dist == "poisson") {
    .do.MCMC.Poisson(fdata, model, prior, mcmc)
  } else if (dist == "cond.poisson") {
    .do.MCMC.CondPoisson(fdata, model, prior, mcmc)
  } else if (dist == "binomial") {
    .do.MCMC.Binomial(fdata, model, prior, mcmc)
  } else if (dist == "exponential") {
    .do.MCMC.Exponential(fdata, model, prior, mcmc)
  } else if (dist == "normal") {
    .do.MCMC.Normal(fdata, model, prior, mcmc)
  } else if (dist == "student") {
    .do.MCMC.Student(fdata, model, prior, mcmc)
  } else if (dist == "normult") {
    .do.MCMC.Normult(fdata, model, prior, mcmc)
  } else if (dist == "studmult") {
    .do.MCMC.Studmult(fdata, model, prior, mcmc)
  }
} ## end mixturemcmc

### Private functions
### These functions are not exported

### Checking
### Check arguments: 'fdata' must contain valid data in @y and in case of
### starting with sampling the parameters indicators in @S. Further,
### the data in @y must match with the specified distribution in @dist
### of 'model'.
### If it should started with sampling the indicators, 'model' must
### contain valid starting parameters in @par and @weight.
### The 'prior' object must contain valid parameters for the prior
### distribution.
### Further, if a fixed indicator model is used, @startpar in 'mcmc'
### must be TRUE and @ranperm must be FALSE.
".check.args.Mixturemcmc" <- function(fdata.obj, model.obj,
                                      prior.obj, mcmc.obj, n.args) {
  ## Check if all arguments are provided
  if (n.args < 4) {
    stop("All arguments must be provided.", call. = FALSE)
  }
  ## Check if 'fdata' object is valid
  if (class(fdata.obj) != "fdata") {
    stop(paste("Unkown argument. Argument 1 must be an ",
      "object of class 'fdata'.",
      sep = ""
    ), call. = FALSE)
  }
  hasY(fdata.obj, verbose = TRUE)
  ## Check if 'model' was provided:
  if (class(model.obj) != "model") {
    stop(paste("Unknown argument. Argument 2 must be an ",
      "object of class 'model'.",
      sep = ""
    ), call. = FALSE)
  }
  ## Check if 'prior' was provided:
  if (class(prior.obj) != "prior") {
    stop(paste("Unknown argument. Argument 3 must be an ",
      "object of class 'prior'.",
      sep = ""
    ), call. = FALSE)
  }
  ## Check if 'mcmc' was provided:
  if (class(mcmc.obj) != "mcmc") {
    stop(paste("Unkown argument. Argument 4 must be an ",
      "object of class 'mcmc'.",
      sep = ""
    ), call. = FALSE)
  }
  ## Check if @startpar in 'mcmc' object and @indicfix in
  ## 'model' object match.
  ## For fixed indicator models indicators are not sampled.
  if (model.obj@indicfix && !mcmc.obj@startpar) {
    mcmc.obj@startpar <- TRUE
  }
  ## Check if @K in 'model' object is one. For a model with
  ## only one component indicators are not sampled.
  if (model.obj@K == 1) {
    mcmc.obj@startpar <- TRUE
  }
  ## If @startpar in 'mcmc.obj' is TRUE, it should be started
  ## by sampling the parameters. In this case starting
  ## indicators must be provided in the 'fdata.obj' object.
  ## If @startpar in 'mcmc.obj' is FALSE it should be started
  ## by sampling the indicators. In this case starting
  ## parameters must be provided in the 'model.obj' object.
  if (model.obj@K > 1) {
    if (mcmc.obj@startpar) {
      if (!hasS(fdata.obj)) {
        stop(paste("For starting with sampling the parameters ",
          "the 'fdata' object must contain starting ",
          "indicator values. See ?mcmcstart for ",
          "generating valid starting values.",
          sep = ""
        ),
        call. = FALSE
        )
      }
    } else {
      if (!hasPar(model.obj)) {
        stop(paste("For starting with sampling the indicators ",
          "the 'model' object must contain starting ",
          "parameter values. See ?mcmcstart for ",
          "generating valid starting values.",
          sep = ""
        ),
        call. = FALSE
        )
      }
      if (!hasWeight(model.obj)) {
        stop(paste("For starting with sampling the indicators ",
          "the 'model' object must contain starting ",
          "weight values. See ?mcmcstart for ",
          "generating valid starting values.",
          sep = ""
        ),
        call. = FALSE
        )
      }
    }
  }
  ## Check if 'fdata' object and 'model' objects match
  ## Call '.check.fdata.model.Mcmcstart()' from 'mcmcstart.R'.
  .check.fdata.model.Mcmcstart(fdata.obj, model.obj)
  ## Check if 'prior' object is valid
  if (!hasPriorPar(prior.obj, model.obj)) {
    stop(paste("Slot @par in 'prior' object is empty. ",
      "For MCMC sampling the prior needs fully ",
      "specified parameters. See ?priordefine for ",
      "generating valid prior parameters.",
      sep = ""
    ),
    call. = FALSE
    )
  }
  if (!model.obj@indicfix && model.obj@K > 1) {
    if (!hasPriorWeight(prior.obj, model.obj)) {
      stop(paste("Slot @weight of 'prior' object is empty. ",
        "For MCMC sampling the prior needs specified ",
        "parameters for the prior of the weights. See ",
        "?priordefine for generating valid prior ",
        "parameters.",
        sep = ""
      ), call. = FALSE)
    }
  }
  ## Check if @indicfix in 'model' object and
  ## @ranperm in 'mcmc' object match.
  ## For a fixed indicator model random permutation
  ## sampling is senseless.
  if (model.obj@indicfix && mcmc.obj@ranperm) {
    mcmc.obj@ranperm <- FALSE
  }
  ## For a model with only one component random permutation
  ## is senseless as well.
  if (model.obj@K == 1 && mcmc.obj@ranperm) {
    mcmc.obj@ranperm <- FALSE
  }
  return(mcmc.obj)
}

### Validity
### For a Binomial model either the 'data' object
### or the 'model' object must have specified
### repetitions 'T'. This can be either a 'matrix'
### object of dimension N x 1 or 1 x 1 (if all
### repetitions are the same)
".valid.Reps.Binomial" <- function(data, model) {
  has.reps <- !all(is.na(data@T))
  if (has.reps) {
    if (data@bycolumn) {
      if (nrow(data@T) != N && nrow(data@T) != 1) {
        stop(paste("Number of repetitions in slot @T of 'data' object ",
          "does not match number of observations in slot @N.",
          sep = ""
        ), call. = FALSE)
      } else if (nrow(data@T) == N) {
        T <- data@T
      } else { ## dimension of T is 1 x 1
        T <- matrix(data@T[1, 1], nrow = N, ncol = 1)
      }
    } else { ## data stored by row
      if (ncol(data@T) != N && ncol(data@T) != 1) {
        stop(paste("Number of repetitions in slot @T of 'data' object ",
          "does not match number of observations slot @N.",
          sep = ""
        ), call. = FALSE)
      } else if (ncol(data@T) == N) {
        T <- t(data@T)
      } else { ## dimension of T is 1 x 1
        T <- matrix(data@T[1, 1], nrow = N, ncol = 1)
      }
    }
  } else { ## then check in model
    has.reps <- !all(is.na(model@T))
    if (has.reps) {
      if (nrow(model@T) != N && nrow(model@T) != 1) {
        stop(paste("Neither 'data' nor 'model' has correctly ",
          "specified repetitions in slot @T for a binomial model.",
          sep = ""
        ), call. = FALSE)
      } else if (nrow(model@T) == N) {
        T <- model@T
      } else { ## dimension of T is 1 x 1
        T <- matrix(model@T[1, 1], nrow = N, ncol = 1)
      }
    } else {
      stop(paste("Neither 'data' object nor 'model' object has ",
        "repetitions in slot @T for a binomial model specified.",
        sep = ""
      ), call. = FALSE)
    }
  }
  ## Check for identifiability ##
  ## Reference: Teicher (1961) ##
  rep.occ <- table(T)
  if (dim(unique(T))[1] == 1) {
    if (T[1, 1] < 2 * model@K - 1) {
      warning(paste("This binomial mixture model is not identifiable. ",
        "For equal repetitions in slot @T it must hold T >= 2K - 1. ",
        "See Teicher (1961) for reference.",
        sep = ""
      ), call. = FALSE)
    }
  } else {
    if (length(rep.occ) != nrow(T)) {
      if (all(dimnames(rep.occ)$T < rep.occ - 1)) {
        warning(paste("This binomial mixture model is not identifiable. ",
          "For varying repetitions 'T_i' in slot @T it must hold T_h ",
          "> r_h - 1, for unique repetitions 'T_h' and their ",
          "respective occurences 'r_h'. See Teicher (1961) ",
          "for reference",
          sep = ""
        ), call. = FALSE)
      } else {
        diff <- diff(sort(unique(T)))
        if (any(diff < rep.occ[1:(length(diff))])) {
          warning(paste("This binomial mixture model is not identifiable. ",
            "For varying repetitions 'T_i' in slot @T it must hold T_h ",
            "- T_(h+1) >= r_h for unique repetitions 'T_h' and ",
            "respective occurrences 'r_h'. See Teicher (1961) ",
            "for reference.",
            sep = ""
          ), call. = FALSE)
        }
      }
    }
  }
}

### MCMC
### For each model the MCMC output has to be prepared
### MCMC Poisson: Prepares all data containers for MCMC sampling for
### Poisson mixture models regarding the specifications in 'prior.obj'
### 'model.obj' and 'mcmc.obj'.
".do.MCMC.Poisson" <- function(fdata.obj, model.obj, prior.obj, mcmc.obj) {
  ## Base slots inherited to every derived class
  K <- model.obj@K
  N <- fdata.obj@N
  M <- mcmc.obj@M
  ranperm <- mcmc.obj@ranperm
  burnin <- mcmc.obj@burnin
  ## Set for MCMC default exposures:
  if (!hasExp(fdata.obj)) {
    fdata.obj@exp <- matrix(1, nrow = N, ncol = 1)
  }
  pars <- list(lambda = array(numeric(), dim = c(M, K)))
  log.mixlik <- array(numeric(), dim = c(M, 1))
  log.mixprior <- array(numeric(), dim = c(M, 1))
  if (mcmc.obj@storepost) {
    post.a <- array(numeric(), dim = c(M, K))
    post.b <- array(numeric(), dim = c(M, K))
    post.par <- list(a = post.a, b = post.b)
    posts <- list(par = post.par)
    if (!model.obj@indicfix) {
      posts$weight <- array(numeric(), dim = c(M, K))
    }
  }
  ## Model with fixed indicators
  if (model.obj@indicfix || K == 1) {
    logs <- list(mixlik = log.mixlik, mixprior = log.mixprior)
    ## Model with simple prior
    if (!prior.obj@hier) {
      ## Model output with NO posterior parameters stored
      if (!mcmc.obj@storepost) {
        mcmcout <- .mcmcoutputfix(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          model = model.obj,
          prior = prior.obj
        )
        .Call("mcmc_poisson_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        return(mcmcout)
      } else {
        ## Model output with posterior parameters stored ##
        mcmcout <- .mcmcoutputfixpost(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          post = posts,
          model = model.obj,
          prior = prior.obj
        )
        .Call("mcmc_poisson_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        return(mcmcout)
      }
      ## end no hier
    } else {
      ## Model with hierarchical prior ##
      hypers <- list(b = array(numeric(), dim = c(M, 1)))
      ## Model output with NO posterior parameters stored ##
      if (!mcmc.obj@storepost) {
        mcmcout <- .mcmcoutputfixhier(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          hyper = hypers,
          model = model.obj,
          prior = prior.obj
        )
        .Call("mcmc_poisson_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        return(mcmcout)
      } else {
        ## Model output with posterior parameters stored ##
        mcmcout <- .mcmcoutputfixhierpost(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          hyper = hypers, post = posts,
          model = model.obj,
          prior = prior.obj
        )
        .Call("mcmc_poisson_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        return(mcmcout)
      }
      ## end hier
    }
    ## end indicfix
  } else if (!model.obj@indicfix && K > 1) {
    ## Model with simulated indicators ##
    log.cdpost <- array(numeric(), dim = c(M, 1))
    logs <- list(
      mixlik = log.mixlik,
      mixprior = log.mixprior,
      cdpost = log.cdpost
    )
    weights <- array(numeric(), dim = c(M, K))
    entropies <- array(numeric(), dim = c(M, 1))
    STm <- array(integer(), dim = c(M, 1))
    Sm <- array(integer(), dim = c(N, mcmc.obj@storeS))
    NKm <- array(integer(), dim = c(M, K))
    clustm <- array(integer(), dim = c(N, 1))
    if (!mcmc.obj@startpar) {
      ## First sample for the indicators
      datac <- dataclass(fdata.obj, model.obj, simS = TRUE)
      Sm[, 1] <- as.integer(datac$S)
    }
    ## Model with simple prior ##
    if (!prior.obj@hier) {
      ## Model output with NO posterior parameters stored ##
      if (!mcmc.obj@storepost) {
        mcmcout <- .mcmcoutputbase(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          weight = weights,
          entropy = entropies,
          ST = STm, S = Sm, NK = NKm,
          clust = clustm, model = model.obj,
          prior = prior.obj
        )
        .Call("mcmc_poisson_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        if (mcmc.obj@storeS == 0) {
          mcmcout@S <- as.array(NA)
        }
        return(mcmcout)
      } else {
        ## Model output with posterior parameters stored ##
        mcmcout <- .mcmcoutputpost(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          weight = weights,
          entropy = entropies,
          ST = STm, S = Sm, NK = NKm,
          clust = clustm, post = posts,
          model = model.obj, prior = prior.obj
        )
        .Call("mcmc_poisson_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        if (mcmc.obj@storeS == 0) {
          mcmcout@S <- as.array(NA)
        }
        return(mcmcout)
      }
      ## end no hier
    } else {
      ## model with hierarchical prior ##
      hypers <- list(b = array(numeric(), dim = c(M, 1)))
      ## model output with NO posterior parameters stored ##
      if (!mcmc.obj@storepost) {
        mcmcout <- .mcmcoutputhier(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          weight = weights,
          entropy = entropies,
          ST = STm, S = Sm, NK = NKm,
          clust = clustm, hyper = hypers,
          model = model.obj, prior = prior.obj
        )
        .Call("mcmc_poisson_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        if (mcmc.obj@storeS == 0) {
          mcmcout@S <- as.array(NA)
        }
        return(mcmcout)
      } else {
        ## model output with posterior parameters stored ##
        mcmcout <- .mcmcoutputhierpost(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          weight = weights,
          entropy = entropies,
          ST = STm, S = Sm, NK = NKm,
          clust = clustm, hyper = hypers,
          post = posts,
          model = model.obj, prior = prior.obj
        )
        .Call("mcmc_poisson_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        if (mcmc.obj@storeS == 0) {
          mcmcout@S <- as.array(NA)
        }
        return(mcmcout)
      }
    } ## end hier
  } ## end no indicfix
}

### ----------------------------------------------------------------------------
### .do.MCMC.Binomial
### @description    Performs MCMC simulation for A Binomial mixture model using
###                 the Gibbs Sampler.
### @par    fdata.obj   an S4 object of class 'fdata'
### @par    model.obj   an S4 object of class 'model'
### @par    prior.obj   an S4 object of class 'prior'
### @par    mcmc.obj    an S4 object of class 'mcmc'
### @return         an S4 object of class 'mcmcoutput'
### @see ?mixturemcmc, ?fdata, ?model, ?prior, ?mcmc, ?mcmcoutput
### @author Lars Simon Zehnder
### ----------------------------------------------------------------------------
".do.MCMC.Binomial" <- function(fdata.obj, model.obj, prior.obj, mcmc.obj) {
  ## Base slots inherited to every derived 'mcmcoutput' class
  K <- model.obj@K
  N <- fdata.obj@N
  M <- mcmc.obj@M
  burnin <- mcmc.obj@burnin
  ranperm <- mcmc.obj@ranperm
  pars <- list(p = array(numeric(), dim = c(c(M, K))))
  log.mixlik <- array(numeric(), dim = c(M, 1))
  log.mixprior <- array(numeric(), dim = c(M, 1))
  if (mcmc.obj@storepost) {
    post.a <- array(numeric(), dim = c(M, K))
    post.b <- array(numeric(), dim = c(M, K))
    post.par <- list(a = post.a, b = post.b)
    posts <- list(par = post.par)
    if (!model.obj@indicfix) {
      posts$weight <- array(numeric(), dim = c(M, K))
    }
  }
  ## Model with fixed indicators
  if (model.obj@indicfix || K == 1) {
    logs <- list(mixlik = log.mixlik, mixprior = log.mixprior)
    if (!mcmc.obj@storepost) {
      mcmcout <- .mcmcoutputfix(
        M = M, burnin = burnin, ranperm = ranperm,
        par = pars, log = logs, model = model.obj,
        prior = prior.obj
      )
      .Call("mcmc_binomial_cc", fdata.obj, model.obj, prior.obj,
        mcmc.obj, mcmcout,
        PACKAGE = "finmix"
      )
      return(mcmcout)
    } else {
      ## MCMC output with posterior hyper parameters stored
      mcmcout <- .mcmcoutputfixpost(
        M = M, burnin = burnin, ranperm = ranperm,
        par = pars, log = logs, post = posts,
        model = model.obj, prior = prior.obj
      )
      .Call("mcmc_binomial_cc", fdata.obj, model.obj, prior.obj,
        mcmc.obj, mcmcout,
        PACKAGE = "finmix"
      )
      return(mcmcout)
    }
    ## End: indicfix
  } else if (!model.obj@indicfix && K > 1) {
    ## Model with simulated indicators
    log.cdpost <- array(numeric(), dim = c(M, 1))
    logs <- list(
      mixlik = log.mixlik, mixprior = log.mixprior,
      cdpost = log.cdpost
    )
    weights <- array(numeric(), dim = c(M, K))
    entropies <- array(numeric(), dim = c(M, 1))
    STm <- array(integer(), dim = c(M, 1))
    Sm <- array(integer(), dim = c(N, mcmc.obj@storeS))
    NKm <- array(integer(), dim = c(M, K))
    clustm <- array(integer(), dim = c(N, 1))
    if (!mcmc.obj@startpar) {
      ## First sample for the indicators
      datac <- dataclass(fdata.obj, model.obj, simS = TRUE)
      Sm[, 1] <- as.integer(datac$S)
    }
    if (!mcmc.obj@storepost) {
      mcmcout <- .mcmcoutputbase(
        M = M, burnin = burnin, ranperm = ranperm,
        par = pars, log = logs, weight = weights,
        entropy = entropies, ST = STm, S = Sm,
        NK = NKm, clust = clustm,
        model = model.obj, prior = prior.obj
      )
      .Call("mcmc_binomial_cc", fdata.obj, model.obj, prior.obj, mcmc.obj,
        mcmcout,
        PACKAGE = "finmix"
      )
      if (mcmc.obj@storeS == 0) {
        mcmcout@S <- as.array(NA)
      }
      return(mcmcout)
    } else {
      ## MCMC output with posterior hyper parameters stored
      mcmcout <- .mcmcoutputpost(
        M = M, burnin = burnin, ranperm = ranperm,
        par = pars, log = logs, weight = weights,
        entropy = entropies, ST = STm, S = Sm,
        NK = NKm, clust = clustm, post = posts,
        model = model.obj, prior = prior.obj
      )
      .Call("mcmc_binomial_cc", fdata.obj, model.obj, prior.obj, mcmc.obj,
        mcmcout,
        PACKAGE = "finmix"
      )
      if (mcmc.obj@storeS == 0) {
        mcmcout@S <- as.array(NA)
      }
      return(mcmcout)
    }
  } ## End no indicfix
}

### -------------------------------------------------------------------------
### .do.MCMC.Exponential
### @description    Prepares all object for the MCMC simulation of an
###                 Exponential model.
### @param  fdata.obj   an S4 object of class 'fdata.obj'
### @param  model.obj   an S4 object of class 'model'
### @param  prior.obj   an S4 object of class 'prior'
### @param  mcmc.obj    an S4 object of class 'mcmc'
### @return an S4 object of class union 'mcmcoutput'
### @detail Internally the C++ routine 'mcmc_exponential_cc' is called
### @see    ?mixturemcmc, mcmc_exponential_cc
### @author Lars Simon Zehnder
### -------------------------------------------------------------------------
".do.MCMC.Exponential" <- function(fdata.obj, model.obj, prior.obj, mcmc.obj) {
  # Base slots inherited to each derived class
  K <- model.obj@K
  N <- fdata.obj@N
  M <- mcmc.obj@M
  ranperm <- FALSE
  burnin <- mcmc.obj@burnin
  pars <- list(lambda = array(numeric(), dim = c(M, K)))
  log.mixlik <- array(numeric(), dim = c(M, 1))
  log.mixprior <- array(numeric(), dim = c(M, 1))
  if (mcmc.obj@storepost) {
    post.a <- array(numeric(), dim = c(M, K))
    post.b <- array(numeric(), dim = c(M, K))
    post.par <- list(a = post.a, b = post.b)
    posts <- list(par = post.par)
    if (!model.obj@indicfix) {
      posts$weight <- array(numeric(), dim = c(M, K))
    }
  }
  # Model with fixed indicators
  if (model.obj@indicfix || K == 1) {
    logs <- list(mixlik = log.mixlik, mixprior = log.mixprior)
    # Model with simple prior
    if (!mcmc.obj@storepost) {
      # Model output with NO posterior parameters stored
      mcmcout <- .mcmcoutputfix(
        M = M, burnin = burnin, ranperm = ranperm,
        par = pars, log = logs, model = model.obj,
        prior = prior.obj
      )
    } else {
      # Model output with posterior parameters stored
      mcmcout <- .mcmcoutputfixpost(
        M = M, burnin = burnin, ranperm = ranperm,
        par = pars, log = logs, post = posts,
        model = model.obj, prior = prior.obj
      )
    }
  } else if (!model.obj@indicfix && K > 1) {
    # Model with simulated indicators
    log.cdpost <- array(numeric(), dim = c(M, 1))
    logs <- list(
      mixlik = log.mixlik, mixprior = log.mixprior,
      cdpost = log.cdpost
    )
    weights <- array(numeric(), dim = c(M, K))
    entropies <- array(numeric(), dim = c(M, 1))
    STm <- array(integer(), dim = c(M, 1))
    Sm <- array(integer(), dim = c(N, mcmc.obj@storeS))
    NKm <- array(integer(), dim = c(M, K))
    clustm <- array(integer(), dim = c(N, 1))
    if (!mcmc.obj@startpar) {
      # First sample for the indicators
      datac <- dataclass(fdata.obj, model.obj, simS = TRUE)
      SM[, 1] <- as.integer(datac$S)
    }
    # Model with simple prior
    # Model output with NO posterior parameters stored
    if (!mcmc.obj@storepost) {
      mcmcout <- .mcmcoutputbase(
        M = M, burnin = burnin, ranperm = ranperm,
        par = pars, log = logs, weight = weights,
        entropy = entropies, ST = STm, S = Sm,
        NK = NKm, clust = clustm, model = model.obj,
        prior = prior.obj
      )
    } else {
      # Model output with posterior parameters stored
      mcmcout <- .mcmcoutputpost(
        M = M, burnin = burnin, ranperm = ranperm,
        par = pars, log = logs, weight = weights,
        entropy = entropies, ST = STm, S = Sm,
        NK = NKm, clust = clustm, post = posts,
        model = model.obj, prior = prior.obj
      )
    }
  }
  .Call("mcmc_exponential_cc", fdata.obj, model.obj, prior.obj, mcmc.obj,
    mcmcout,
    PACKAGE = "finmix"
  )
  if (mcmc.obj@storeS == 0) {
    mcmcout@S <- as.array(NA)
  }
  return(mcmcout)
}

".do.MCMC.CondPoisson" <- function(fdata.obj, model.obj, prior.obj, mcmc.obj) {
  if (nrow(fdata.obj@exp) == 1) {
    if (is.na(fdata.obj@exp)) {
      fdata.obj@exp <- matrix(1, nrow = fdata.obj@N, ncol = 1)
    } else {
      fdata.obj@exp <- matrix(fdata.obj@exp, nrow = fdata.obj@N, ncol = 1)
    }
  }
  if (mcmc.obj@ranperm) {
    mcmc.obj@ranperm <- FALSE
  }
  ## base slots inherited to every derived class ##
  K <- model.obj@K
  N <- fdata.obj@N
  M <- mcmc.obj@M
  burnin <- mcmc.obj@burnin
  ranperm <- mcmc.obj@ranperm
  pars <- list(
    lambda = array(numeric(), dim = c(M, K)),
    acc = 0.0
  )
  log.mixlik <- array(numeric(), dim = c(M, 1))
  log.mixprior <- array(numeric(), dim = c(M, 1))
  if (mcmc.obj@storepost) {
    post.Q <- array(numeric(), dim = c(M, K))
    post.N <- array(numeric(), dim = c(M, K))
    post.par <- list(Q = post.Q, N = post.N)
    posts <- list(par = post.par)
    if (!model.obj@indicfix) {
      posts$weight <- array(numeric(), dim = c(M, K))
    }
  }
  ## model with fixed indicators ##
  if (model.obj@indicfix || K == 1) {
    logs <- list(mixlik = log.mixlik, mixprior = log.mixprior)
    ## model with simple prior ##
    if (!prior.obj@hier) {
      ## model output with NO posterior parameters stored ##
      if (!mcmc.obj@storepost) {
        mcmcout <- .mcmcoutputfix(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          model = model.obj,
          prior = prior.obj
        )
        .Call("mcmc_condpoisson_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        return(mcmcout)
      } else {
        ## model output with posterior parameters stored ##
        mcmcout <- .mcmcoutputfixpost(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          post = posts,
          model = model.obj,
          prior = prior.obj
        )
        .Call("mcmc_condpoisson_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        return(mcmcout)
      }
      ## end no hier
    } else {
      ## model with hierarchical prior ##
      hypers <- list(b = array(numeric(), dim = c(M, 1)))
      ## model output with NO posterior parameters stored ##
      if (!mcmc.obj@storepost) {
        mcmcout <- .mcmcoutputfixhier(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          hyper = hypers,
          model = model.obj,
          prior = prior.obj
        )
        .Call("mcmc_condpoisson_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        return(mcmcout)
      } else {
        ## model output with posterior parameters stored ##
        mcmcout <- .mcmcoutputfixhierpost(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          hyper = hypers,
          post = posts,
          model = model.obj,
          prior = prior.obj
        )
        .Call("mcmc_condpoisson_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        return(mcmcout)
      }
      ## end hier
    }
    ## end indicfix
  } else if (!model.obj@indicfix && K > 1) {
    ## model with simulated indicators ##
    log.cdpost <- array(numeric(), dim = c(M, 1))
    logs <- list(
      mixlik = log.mixlik, mixprior = log.mixprior,
      cdpost = log.cdpost
    )
    weights <- array(numeric(), dim = c(M, K))
    entropies <- array(numeric(), dim = c(M, 1))
    STm <- array(integer(), dim = c(M, 1))
    Sm <- array(integer(), dim = c(N, mcmc.obj@storeS))
    NKm <- array(integer(), dim = c(M, K))
    clustm <- array(integer(), dim = c(N, 1))
    if (!mcmc.obj@startpar) {
      ## First sample for the indicators
      datac <- dataclass(fdata.obj, model.obj, simS = TRUE)
      Sm[, 1] <- as.integer(datac$S)
    }
    ## model with simple prior ##
    if (!prior.obj@hier) {
      ## model output with NO posterior parameters stored ##
      if (!mcmc.obj@storepost) {
        mcmcout <- .mcmcoutputbase(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          weight = weights,
          entropy = entropies,
          ST = STm, S = Sm, NK = NKm,
          clust = clustm,
          model = model.obj, prior = prior.obj
        )
        .Call("mcmc_condpoisson_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        if (mcmc.obj@storeS == 0) {
          mcmcout@S <- as.array(NA)
        }
        return(mcmcout)
      } else {
        ## model output with posterior parameters stored ##
        mcmcout <- .mcmcoutputpost(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          weight = weights,
          entropy = entropies,
          ST = STm, S = Sm, NK = NKm,
          clust = clustm, post = posts,
          model = model.obj, prior = prior.obj
        )
        .Call("mcmc_condpoisson_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        if (mcmc.obj@storeS == 0) {
          mcmcout@S <- as.array(NA)
        }
        return(mcmcout)
      }
      ## end no hier
    } else {
      ## model with hierarchical prior ##
      hypers <- list(b = array(numeric(), dim = c(M, 1)))
      ## model output with NO posterior parameters stored ##
      if (!mcmc.obj@storepost) {
        mcmcout <- .mcmcoutputhier(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          weight = weights,
          entropy = entropies,
          ST = STm, S = Sm, NK = NKm,
          clust = clustm, hyper = hypers,
          model = model.obj, prior = prior.obj
        )
        .Call("mcmc_condpoisson_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        if (mcmc.obj@storeS == 0) {
          mcmcout@S <- as.array(NA)
        }
        return(mcmcout)
      } else {
        ## model output with posterior parameters stored ##
        mcmcout <- .mcmcoutputhierpost(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          weight = weights,
          entropy = entropies,
          ST = STm, S = Sm, NK = NKm,
          clust = clustm, hyper = hypers,
          post = posts,
          model = model.obj, prior = prior.obj
        )
        .Call("mcmc_condpoisson_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        if (mcmc.obj@storeS == 0) {
          mcmcout@S <- as.array(NA)
        }
        return(mcmcout)
      }
    } ## end hier
  } ## end no indicfix
}

".do.MCMC.Normal" <- function(fdata.obj, model.obj, prior.obj,
                              mcmc.obj) {
  ## Base slots inherited to each derived class
  K <- model.obj@K
  N <- fdata.obj@N
  M <- mcmc.obj@M
  ranperm <- mcmc.obj@ranperm
  burnin <- mcmc.obj@burnin
  ## Set for MCMC default exposures:
  pars <- list(
    mu = array(numeric(), dim = c(M, K)),
    sigma = array(numeric(), dim = c(M, K))
  )
  log.mixlik <- array(numeric(), dim = c(M, 1))
  log.mixprior <- array(numeric(), dim = c(M, 1))
  if (mcmc.obj@storepost) {
    b.post <- array(numeric(), dim = c(M, K))
    B.post <- array(numeric(), dim = c(M, K))
    mu.post <- list(b = b.post, B = B.post)
    c.post <- array(numeric(), dim = c(M, K))
    C.post <- array(numeric(), dim = c(M, K))
    sigma.post <- list(c = c.post, C = C.post)
    par.post <- list(mu = mu.post, sigma = sigma.post)
    posts <- list(par = par.post)
    if (!model.obj@indicfix) {
      posts$weight <- array(numeric(), dim = c(M, K))
    }
  }
  ## Model with fixed indicators
  if (model.obj@indicfix || K == 1) {
    logs <- list(mixlik = log.mixlik, mixprior = log.mixprior)
    ## Model with simple prior
    if (!prior.obj@hier) {
      ## Model output with NO posterior parameters stored
      if (!mcmc.obj@storepost) {
        mcmcout <- .mcmcoutputfix(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          model = model.obj,
          prior = prior.obj
        )
        .Call("mcmc_normal_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        return(mcmcout)
      } else {
        mcmcout <- .mcmcoutputfixpost(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          post = posts,
          model = model.obj,
          prior = prior.obj
        )
        .Call("mcmc_normal_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        return(mcmcout)
      }
      ## end no hier
    } else {
      ## Model with hierarchical prior
      hypers <- list(C = array(numeric(), dim = c(M, 1)))
      ## Model with NO posterior parameters stored
      if (!mcmc.obj@storepost) {
        mcmcout <- mcmcoutputfixhier(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          hyper = hypers,
          model = model.obj,
          prior = prior.obj
        )
        .Call("mcmc_normal_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        return(mcmcout)
      } else {
        ## Model with posterior parameters stored
        mcmcout <- mcmcoutputfixhierpost(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logd,
          hyper = hypers, post = posts,
          model = model.obj,
          prior = prior.obj
        )
        .Call("mcmc_normal_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        return(mcmcout)
      }
      ## end hier
    }
    ## end indicfix
  } else if (!model.obj@indicfix && K > 1) {
    ## Model with simulated indicators
    log.cdpost <- array(numeric(), dim = c(M, 1))
    logs <- list(
      mixlik = log.mixlik,
      mixprior = log.mixprior,
      cdpost = log.cdpost
    )
    weights <- array(numeric(), dim = c(M, K))
    entropies <- array(numeric(), dim = c(M, 1))
    STm <- array(integer(), dim = c(M, 1))
    Sm <- array(integer(), dim = c(N, mcmc.obj@storeS))
    NKm <- array(integer(), dim = c(M, K))
    clustm <- array(integer(), dim = c(N, 1))
    if (!mcmc.obj@startpar) {
      ## First sample for the indicators
      datac <- dataclass(fdata.obj, model.obj, simS = TRUE)
      Sm[, 1] <- as.integer(datac$s)
    }
    ## Model with simple prior
    if (!prior.obj@hier) {
      ## Model output with NO posterior parameters stored
      if (!mcmc.obj@storepost) {
        mcmcout <- .mcmcoutputbase(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          weight = weights,
          entropy = entropies,
          ST = STm, S = Sm, NK = NKm,
          clust = clustm, model = model.obj,
          prior = prior.obj
        )
        .Call("mcmc_normal_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        if (mcmc.obj@storeS == 0) {
          mcmcout@S <- as.array(as.integer(NA))
        }
        return(mcmcout)
      } else {
        ## Model output with posterior parameters stored
        mcmcout <- .mcmcoutputpost(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          weight = weights,
          entropy = entropies,
          ST = STm, S = Sm, NK = NKm,
          clust = clustm, post = posts,
          model = model.obj, prior = prior.obj
        )
        .Call("mcmc_normal_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        if (mcmc.obj@storeS == 0) {
          mcmcout@S <- as.array(as.integer(NA))
        }
        return(mcmcout)
      }
      ## end no hier
    } else {
      ## Model with hierarchical prior
      hypers <- list(C = array(numeric(), dim = c(M, 1)))
      ## Model with NO posterior parameters stored
      if (!mcmc.obj@storepost) {
        mcmcout <- .mcmcoutputhier(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          weight = weights,
          entropy = entropies,
          ST = STm, S = Sm, NK = NKm,
          clust = clustm, hyper = hypers,
          model = model.obj, prior = prior.obj
        )
        .Call("mcmc_normal_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        if (mcmc.obj@storeS == 0) {
          mcmcout@S <- as.array(ias.integer(NA))
        }
        return(mcmcout)
      } else {
        ## Model with posterior parameters stored
        mcmcout <- .mcmcoutputhierpost(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          weight = weights,
          entropy = entropies,
          ST = STm, S = Sm, NK = NKm,
          clust = clustm, hyper = hypers,
          post = posts, model = model.obj,
          prior = prior.obj
        )
        .Call("mcmc_normal_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        if (mcmc.obj@storeS == 0) {
          mcmcout@S <- as.array(as.integer(NA))
        }
        return(mcmcout)
      }
    } ## end hier
  } ## end no indicfix
}

".do.MCMC.Student" <- function(fdata.obj, model.obj, prior.obj,
                               mcmc.obj) {
  ## Base slots inherited to each derived class
  K <- model.obj@K
  N <- fdata.obj@N
  M <- mcmc.obj@M
  ranperm <- mcmc.obj@ranperm
  burnin <- mcmc.obj@burnin
  ## Set for MCMC default exposures:
  pars <- list(
    mu = array(numeric(), dim = c(M, K)),
    sigma = array(numeric(), dim = c(M, K)),
    df = array(numeric(), dim = c(M, K)),
    acc = array(0.0, dim = c(1, K))
  )
  log.mixlik <- array(numeric(), dim = c(M, 1))
  log.mixprior <- array(numeric(), dim = c(M, 1))
  if (mcmc.obj@storepost) {
    b.post <- array(numeric(), dim = c(M, K))
    B.post <- array(numeric(), dim = c(M, K))
    mu.post <- list(b = b.post, B = B.post)
    c.post <- array(numeric(), dim = c(M, K))
    C.post <- array(numeric(), dim = c(M, K))
    sigma.post <- list(c = c.post, C = C.post)
    par.post <- list(mu = mu.post, sigma = sigma.post)
    posts <- list(par = par.post)
    if (!model.obj@indicfix) {
      posts$weight <- array(numeric(), dim = c(M, K))
    }
  }
  ## Model with fixed indicators
  if (model.obj@indicfix || K == 1) {
    logs <- list(mixlik = log.mixlik, mixprior = log.mixprior)
    ## Model with simple prior
    if (!prior.obj@hier) {
      ## Model output with NO posterior parameters stored
      if (!mcmc.obj@storepost) {
        mcmcout <- .mcmcoutputfix(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          model = model.obj,
          prior = prior.obj
        )
        .Call("mcmc_student_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        return(mcmcout)
      } else {
        mcmcout <- .mcmcoutputfixpost(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          post = posts,
          model = model.obj,
          prior = prior.obj
        )
        .Call("mcmc_student_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        return(mcmcout)
      }
      ## end no hier
    } else {
      ## Model with hierarchical prior
      hypers <- list(C = array(numeric(), dim = c(M, 1)))
      ## Model with NO posterior parameters stored
      if (!mcmc.obj@storepost) {
        mcmcout <- mcmcoutputfixhier(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          hyper = hypers,
          model = model.obj,
          prior = prior.obj
        )
        .Call("mcmc_student_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        return(mcmcout)
      } else {
        ## Model with posterior parameters stored
        mcmcout <- mcmcoutputfixhierpost(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logd,
          hyper = hypers, post = posts,
          model = model.obj,
          prior = prior.obj
        )
        .Call("mcmc_student_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        return(mcmcout)
      }
      ## end hier
    }
    ## end indicfix
  } else if (!model.obj@indicfix && K > 1) {
    ## Model with simulated indicators
    log.cdpost <- array(numeric(), dim = c(M, 1))
    logs <- list(
      mixlik = log.mixlik,
      mixprior = log.mixprior,
      cdpost = log.cdpost
    )
    weights <- array(numeric(), dim = c(M, K))
    entropies <- array(numeric(), dim = c(M, 1))
    STm <- array(integer(), dim = c(M, 1))
    Sm <- array(integer(), dim = c(N, mcmc.obj@storeS))
    NKm <- array(integer(), dim = c(M, K))
    clustm <- array(integer(), dim = c(N, 1))
    if (!mcmc.obj@startpar) {
      ## First sample for the indicators
      datac <- dataclass(fdata.obj, model.obj, simS = TRUE)
      Sm[, 1] <- as.integer(datac$s)
    }
    ## Model with simple prior
    if (!prior.obj@hier) {
      ## Model output with NO posterior parameters stored
      if (!mcmc.obj@storepost) {
        mcmcout <- .mcmcoutputbase(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          weight = weights,
          entropy = entropies,
          ST = STm, S = Sm, NK = NKm,
          clust = clustm, model = model.obj,
          prior = prior.obj
        )
        .Call("mcmc_student_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        if (mcmc.obj@storeS == 0) {
          mcmcout@S <- as.array(as.integer(NA))
        }
        return(mcmcout)
      } else {
        ## Model output with posterior parameters stored
        mcmcout <- .mcmcoutputpost(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          weight = weights,
          entropy = entropies,
          ST = STm, S = Sm, NK = NKm,
          clust = clustm, post = posts,
          model = model.obj, prior = prior.obj
        )
        .Call("mcmc_student_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        if (mcmc.obj@storeS == 0) {
          mcmcout@S <- as.array(as.integer(NA))
        }
        return(mcmcout)
      }
      ## end no hier
    } else {
      ## Model with hierarchical prior
      hypers <- list(C = array(numeric(), dim = c(M, 1)))
      ## Model with NO posterior parameters stored
      if (!mcmc.obj@storepost) {
        mcmcout <- .mcmcoutputhier(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          weight = weights,
          entropy = entropies,
          ST = STm, S = Sm, NK = NKm,
          clust = clustm, hyper = hypers,
          model = model.obj, prior = prior.obj
        )
        .Call("mcmc_student_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        if (mcmc.obj@storeS == 0) {
          mcmcout@S <- as.array(as.integer(NA))
        }
        return(mcmcout)
      } else {
        ## Model with posterior parameters stored
        mcmcout <- .mcmcoutputhierpost(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          weight = weights,
          entropy = entropies,
          ST = STm, S = Sm, NK = NKm,
          clust = clustm, hyper = hypers,
          post = posts, model = model.obj,
          prior = prior.obj
        )
        .Call("mcmc_student_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        if (mcmc.obj@storeS == 0) {
          mcmcout@S <- as.array(as.integer(NA))
        }
        return(mcmcout)
      }
    } ## end hier
  } ## end no indicfix
}

".do.MCMC.Normult" <- function(fdata.obj, model.obj, prior.obj, mcmc.obj) {
  ## Base slots inherited to each derived class
  K <- model.obj@K
  N <- fdata.obj@N
  M <- mcmc.obj@M
  ranperm <- mcmc.obj@ranperm
  burnin <- mcmc.obj@burnin

  ## Constants simplifying construction
  r <- fdata.obj@r
  s <- r * (r + 1) / 2

  ## Set for MCMC default expousres
  pars <- list(
    mu = array(numeric(), dim = c(M, r, K)),
    sigma = array(numeric(), dim = c(M, s, K)),
    storeinv = mcmc.obj@storeinv
  )
  if (mcmc.obj@storeinv) {
    pars$sigmainv <- array(numeric(), dim = c(M, s, K))
  }
  log.mixlik <- array(numeric(), dim = c(M, 1))
  log.mixprior <- array(numeric(), dim = c(M, 1))
  if (mcmc.obj@storepost) {
    b.post <- array(numeric(), dim = c(M, r, K))
    B.post <- array(numeric(), dim = c(M, s, K))
    mu.post <- list(b = b.post, B = B.post)
    c.post <- array(numeric(), dim = c(M, K))
    C.post <- array(numeric(), dim = c(M, s, K))
    sigma.post <- list(c = c.post, C = C.post)
    par.post <- list(mu = mu.post, sigma = sigma.post)
    posts <- list(par = par.post)
    if (!model.obj@indicfix) {
      posts$weight <- array(numeric(), dim = c(M, K))
    }
  }
  ## Model with fixed indicators
  if (model.obj@indicfix || K == 1) {
    logs <- list(mixlik = log.mixlik, mixprior = log.mixprior)
    ## Model with simple prior
    if (!prior.obj@hier) {
      ## Model output with NO posterior parameters stored
      if (!mcmc.obj@storepost) {
        mcmcout <- .mcmcoutputfix(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          model = model.obj,
          prior = prior.obj
        )
        .Call("mcmc_normult_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        mcmcout@par$storeinv <- NULL
        return(mcmcout)
      } else {
        mcmcout <- .mcmcoutputfixpost(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          post = posts,
          model = model.obj,
          prior = prior.obj
        )
        .Call("mcmc_normult_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        mcmcout@par$storeinv <- NULL
        return(mcmcout)
      }
      ## end no hier
    } else {
      ## Model with hierarchical prior
      hypers <- list(C = array(numeric(), dim = c(M, s)))
      ## Model with NO posterior parameters stored
      if (!mcmc.obj@storepost) {
        mcmcout <- .mcmcoutputfixhier(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          hyper = hypers,
          model = model.obj,
          prior = prior.obj
        )
        .Call("mcmc_normult_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        mcmcout@par$storeinv <- NULL
        return(mcmcout)
      } else {
        ## Model with posterior parameters stored
        mcmcout <- .mcmcoutputfixhierpost(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          hyper = hypers, post = posts,
          model = model.obj,
          prior = prior.obj
        )
        .Call("mcmc_normult_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        mcmcout@par$storeinv <- NULL
        return(mcmcout)
      }
      ## end hier
    }
    ## end indicfix
  } else if (!model.obj@indicfix && K > 1) {
    ## Model with simulated indicators
    log.cdpost <- array(numeric(), dim = c(M, 1))
    logs <- list(
      mixlik = log.mixlik,
      mixprior = log.mixprior,
      cdpost = log.cdpost
    )
    weights <- array(numeric(), dim = c(M, K))
    entropies <- array(numeric(), dim = c(M, 1))
    STm <- array(integer(), dim = c(M, 1))
    Sm <- array(integer(), dim = c(N, mcmc.obj@storeS))
    NKm <- array(integer(), dim = c(M, K))
    clustm <- array(integer(), dim = c(N, 1))
    if (!mcmc.obj@startpar) {
      ## First sample for the indicators
      datac <- dataclass(fdata.obj, model.obj, simS = TRUE)
      Sm[, 1] <- as.integer(datac$s)
    }
    ## Model with simple prior
    if (!prior.obj@hier) {
      ## Model with NO posterior parameters stored
      if (!mcmc.obj@storepost) {
        mcmcout <- .mcmcoutputbase(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          weight = weights,
          entropy = entropies,
          ST = STm, S = Sm, NK = NKm,
          clust = clustm, model = model.obj,
          prior = prior.obj
        )
        .Call("mcmc_normult_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        if (mcmc.obj@storeS == 0) {
          mcmcout@S <- as.array(as.integer(NA))
        }
        mcmcout@par$storeinv <- NULL
        return(mcmcout)
      } else {
        ## Model output with posterior parameters stored
        mcmcout <- .mcmcoutpupost(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          weight = weights,
          entropy = entropies,
          ST = STm, S = Sm, NK = NKm,
          clust = clustm, post = postm,
          model = model.obj, prior = prior.obj
        )
        .Call("mcmc_normult_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        if (mcmc.obj@storeS == 0) {
          mcmcout@S <- as.array(as.integer(NA))
        }
        mcmcout@par$storeinv <- NULL
        return(mcmcout)
      }
      ## end no hier
    } else {
      ## Model with hierarchical prior
      hypers <- list(C = array(numeric(), dim = c(M, s)))
      ## Model with NO posterior parameters stored
      if (!mcmc.obj@storepost) {
        mcmcout <- .mcmcoutputhier(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          weight = weights,
          entropy = entropies,
          ST = STm, S = Sm, NK = NKm,
          clust = clustm, hyper = hypers,
          model = model.obj, prior = prior.obj
        )
        .Call("mcmc_normult_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        if (mcmc.obj@storeS == 0) {
          mcmcout@S <- as.array(as.integer(NA))
        }
        mcmcout@par$storeinv <- NULL
        return(mcmcout)
      } else {
        ## Model with posterior parameters stored
        mcmcout <- .mcmcoutputhierpost(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          weight = weights,
          entropy = entropies,
          ST = STm, S = Sm, NK = NKm,
          clust = clustm, hyper = hypers,
          post = posts, model = model.obj,
          prior = prior.obj
        )
        .Call("mcmc_normult_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        if (mcmc.obj@storeS == 0) {
          mcmcout@S <- as.array(as.integer(NA))
        }
        mcmcout@par$storeinv <- NULL
        return(mcmcout)
      }
    } ## end hier
  } ## end no indicfix
}

".do.MCMC.Studmult" <- function(fdata.obj, model.obj, prior.obj, mcmc.obj) {
  ## Base slots inherited to each derived class
  K <- model.obj@K
  N <- fdata.obj@N
  M <- mcmc.obj@M
  ranperm <- mcmc.obj@ranperm
  burnin <- mcmc.obj@burnin

  ## Constants simplifying construction
  r <- fdata.obj@r
  s <- r * (r + 1) / 2

  ## Set for MCMC default expousres
  pars <- list(
    mu = array(numeric(), dim = c(M, r, K)),
    sigma = array(numeric(), dim = c(M, s, K)),
    df = array(numeric(), dim = c(M, K)),
    acc = array(0.0, dim = c(1, K)),
    storeinv = mcmc.obj@storeinv
  )
  if (mcmc.obj@storeinv) {
    pars$sigmainv <- array(numeric(), dim = c(M, s, K))
  }
  log.mixlik <- array(numeric(), dim = c(M, 1))
  log.mixprior <- array(numeric(), dim = c(M, 1))
  if (mcmc.obj@storepost) {
    b.post <- array(numeric(), dim = c(M, r, K))
    B.post <- array(numeric(), dim = c(M, s, K))
    mu.post <- list(b = b.post, B = B.post)
    c.post <- array(numeric(), dim = c(M, K))
    C.post <- array(numeric(), dim = c(M, s, K))
    sigma.post <- list(c = c.post, C = C.post)
    par.post <- list(mu = mu.post, sigma = sigma.post)
    posts <- list(par = par.post)
    if (!model.obj@indicfix) {
      posts$weight <- array(numeric(), dim = c(M, K))
    }
  }
  ## Model with fixed indicators
  if (model.obj@indicfix || K == 1) {
    logs <- list(mixlik = log.mixlik, mixprior = log.mixprior)
    ## Model with simple prior
    if (!prior.obj@hier) {
      ## Model output with NO posterior parameters stored
      if (!mcmc.obj@storepost) {
        mcmcout <- .mcmcoutputfix(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          model = model.obj,
          prior = prior.obj
        )
        .Call("mcmc_studmult_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        mcmcout@par$storeinv <- NULL
        return(mcmcout)
      } else {
        mcmcout <- .mcmcoutputfixpost(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          post = posts,
          model = model.obj,
          prior = prior.obj
        )
        .Call("mcmc_studmult_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        mcmcout@par$storeinv <- NULL
        return(mcmcout)
      }
      ## end no hier
    } else {
      ## Model with hierarchical prior
      hypers <- list(C = array(numeric(), dim = c(M, s)))
      ## Model with NO posterior parameters stored
      if (!mcmc.obj@storepost) {
        mcmcout <- .mcmcoutputfixhier(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          hyper = hypers,
          model = model.obj,
          prior = prior.obj
        )
        .Call("mcmc_studmult_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        mcmcout@par$storeinv <- NULL
        return(mcmcout)
      } else {
        ## Model with posterior parameters stored
        mcmcout <- .mcmcoutputfixhierpost(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          hyper = hypers, post = posts,
          model = model.obj,
          prior = prior.obj
        )
        .Call("mcmc_studmult_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        mcmcout@par$storeinv <- NULL
        return(mcmcout)
      }
      ## end hier
    }
    ## end indicfix
  } else if (!model.obj@indicfix && K > 1) {
    ## Model with simulated indicators
    log.cdpost <- array(numeric(), dim = c(M, 1))
    logs <- list(
      mixlik = log.mixlik,
      mixprior = log.mixprior,
      cdpost = log.cdpost
    )
    weights <- array(numeric(), dim = c(M, K))
    entropies <- array(numeric(), dim = c(M, 1))
    STm <- array(integer(), dim = c(M, 1))
    Sm <- array(integer(), dim = c(N, mcmc.obj@storeS))
    NKm <- array(integer(), dim = c(M, K))
    clustm <- array(integer(), dim = c(N, 1))
    if (!mcmc.obj@startpar) {
      ## First sample for the indicators
      datac <- dataclass(fdata.obj, model.obj, simS = TRUE)
      Sm[, 1] <- as.integer(datac$s)
    }
    ## Model with simple prior
    if (!prior.obj@hier) {
      ## Model with NO posterior parameters stored
      if (!mcmc.obj@storepost) {
        mcmcout <- .mcmcoutputbase(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          weight = weights,
          entropy = entropies,
          ST = STm, S = Sm, NK = NKm,
          clust = clustm, model = model.obj,
          prior = prior.obj
        )
        .Call("mcmc_studmult_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        if (mcmc.obj@storeS == 0) {
          mcmcout@S <- as.array(as.integer(NA))
        }
        mcmcout@par$storeinv <- NULL
        return(mcmcout)
      } else {
        ## Model output with posterior parameters stored
        mcmcout <- .mcmcoutpupost(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          weight = weights,
          entropy = entropies,
          ST = STm, S = Sm, NK = NKm,
          clust = clustm, post = postm,
          model = model.obj, prior = prior.obj
        )
        .Call("mcmc_studmult_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        if (mcmc.obj@storeS == 0) {
          mcmcout@S <- as.array(as.integer(NA))
        }
        mcmcout@par$storeinv <- NULL
        return(mcmcout)
      }
      ## end no hier
    } else {
      ## Model with hierarchical prior
      hypers <- list(C = array(numeric(), dim = c(M, s)))
      ## Model with NO posterior parameters stored
      if (!mcmc.obj@storepost) {
        mcmcout <- .mcmcoutputhier(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          weight = weights,
          entropy = entropies,
          ST = STm, S = Sm, NK = NKm,
          clust = clustm, hyper = hypers,
          model = model.obj, prior = prior.obj
        )
        .Call("mcmc_studmult_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        if (mcmc.obj@storeS == 0) {
          mcmcout@S <- as.array(as.integer(NA))
        }
        mcmcout@par$storeinv <- NULL
        return(mcmcout)
      } else {
        ## Model with posterior parameters stored
        mcmcout <- .mcmcoutputhierpost(
          M = M, burnin = burnin,
          ranperm = ranperm,
          par = pars, log = logs,
          weight = weights,
          entropy = entropies,
          ST = STm, S = Sm, NK = NKm,
          clust = clustm, hyper = hypers,
          post = posts, model = model.obj,
          prior = prior.obj
        )
        .Call("mcmc_studmult_cc", fdata.obj, model.obj, prior.obj,
          mcmc.obj, mcmcout,
          PACKAGE = "finmix"
        )
        if (mcmc.obj@storeS == 0) {
          mcmcout@S <- as.array(as.integer(NA))
        }
        mcmcout@par$storeinv <- NULL
        return(mcmcout)
      }
    } ## end hier
  } ## end no indicfix
}
