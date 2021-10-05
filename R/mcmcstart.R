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

#' Finmix default starting values
#' 
#' @description 
#' Calling [mcmcstart()] creates starting values for MCMC sampling. Starting 
#' values are constructed for the indicators in the `fdata` argument and the 
#' parameters in the `model` argument. In addition an `mcmc` object can be 
#' provided in the `varargin` argument to set up all slots consistently for a 
#' non-default setting of hyper-parameters.
#' 
#' To assing the returned objects directly to existing names the assignment 
#' operator `%%=%%` can be used together with a formula concatenating each name 
#' with a tilde `~`. See the examples.
#' 
#' @param fdata An `fdata` object containing the data.
#' @param model A `model` object specifying the finite mixture model to be
#'   estimated. 
#' @param varargin Either `NULL` or an `mcmc` object defining (possibly 
#'   non-default) hyper-parameters. If not provided a default `mcmc` object is 
#'   created internally and returned.
#' @return A list containing the `fdata` object, the `model` object and an 
#'   `mcmc` object all set up for MCMC sampling.
#' @export
#' @name mcmcstart
#' 
#' @examples 
#' # Specify a Poisson model.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 0.7)), K = 2)
#' # Simulate data from the model.
#' f_data <- simulate(f_model)
#' # Set up all objects for MCMC sampling.
#' (f_data ~ f_model ~ f_mcmc) %=% mcmcstart(f_data, f_model)
#' 
#' @seealso 
#' * [fdata][fdata_class] for the definition of the `fdata` class
#' * [model][model_class] for the definition of the `model` class
#' * [mcmc][mcmc_class] for the definition of the `mcmc` class
#' * [mixturemcmc()] for the starting MCMC sampling
"mcmcstart" <- function(fdata, model, varargin) {
  ## Check arguments
  .check.fdata.model.Mcmcstart(fdata, model)
  ## Check if mcmc object was given in arguments
  if (nargs() == 2) {
    mcmc <- mcmc()
  } else {
    .check.mcmc.Mcmcstart(varargin)
    mcmc <- varargin
  }
  K <- model@K
  dist <- model@dist
  ## If @startpar is
  ## TRUE (default):      Start by sampling the parameters
  ##                          -> it needs starting indicators S
  ##                      In this case staring indicators are generated
  ##                      by kmeans-clustering independently of
  ##                      the model.
  ## FALSE:               Start by sampling the indicators
  ##                          -> it needs starting parameters model@par
  ##                      In this case starting parameters are generated
  ##                      dependent on the data in @y of 'fdata' and the
  ##                      the model.
  ## If the model has fixed indicators (@indicfix = TRUE), no indicators
  ## are generated.
  if (!model@indicfix) {
    if (mcmc@startpar) {
      if (K > 1) {
        fdata <- .indicators.Mcmcstart(fdata, model)
      }
    } else {
      model <- .parameters.Mcmcstart(fdata, model, mcmc)
    }
  } else {
    warning(paste("Slot 'indicfix' of 'model' object is ",
      "set to TRUE. 'mcmcstart()' does not ",
      "generate indicators nor starting parameters ",
      "for models with fixed indicators.",
      sep = ""
    ))
  }
  obj.list <- list(fdata = fdata, model = model, mcmc = mcmc)
  return(obj.list)
}

### Private functions.
### These functions are not exported.

### Checking
### Check fdata/model: 'fdata' must be a valid 'fdata' object and 'model'
### must be a valid 'model' object. Furthermore, 'fdata' must have a
### non-empty data slot @y.
### If the distributions in 'model' do not correspond to the dimensions
### @r in 'fdata' an error is thrown.
#' Check argument in `mcmcstart`
#' 
#' For internal usage only. This function checks the input arguments `fdata` 
#' and `model` for consistency. Consistency has to be ensured for the slots 
#' `@@dist` in the `model` object and the slot `@@r` in the `fdata` object. 
#' A dimension `r>1` calls for a multivariate distribution specified in the 
#' `model` object. Furthermore, the `fdata` object must contain data in its 
#' `@@y` slot.
#' 
#' @param fdata_obj An `fdata` object containing the data.
#' @param model_obj A `model` object specifying the finite mixture model.
#' @return None. If any check is not passed an error is thrown.
#' @noRd
#' @keywords internal
#' 
#' @seealso 
#' * [fdata][fdata_class] for the definition of the `fdata` class
#' * [model][model_class] for the definition of the `model` class
#' * [mcmc][mcmc_class] for the definition of the `mcmc` class
#' * [mixturemcmc()] for the starting MCMC sampling
".check.fdata.model.Mcmcstart" <- function(fdata.obj, model.obj) {
  .valid.Fdata(fdata.obj)
  .valid.Model(model.obj)
  hasY(fdata.obj, verbose = TRUE)
  if (fdata.obj@r > 1 && model.obj@dist %in% .get.univ.Model()) {
    stop(paste("Wrong specification of slot 'r' in 'fdata' object. ",
      "Univariate distribution in slot 'dist' of 'model' ",
      "object but dimension in slot 'r' of 'fdata' object ",
      "greater 1.",
      sep = ""
    ))
  } else if (fdata.obj@r < 2 && model.obj@dist %in% .get.multiv.Model()) {
    stop(paste("Wrong specification of slot 'r' ind 'fdata' object ",
      "Multivariate distribution in slot 'dist' if 'model' ",
      "object but dimension in slot 'r' of 'fdata' object ",
      "less than two.",
      sep = ""
    ))
  }
}

### Check varargin: Argument 'varargin' must be an object of class 'mcmc'.
#' Check argument `varargin` in `mcmcstart`
#' 
#' For internal usage only. This function checks the `varargin` input argument. 
#' More specifically, it checks if the argument is an `mcmc` object and if it 
#' is correctly specified.
#' 
#' @param fdata_obj An `fdata` object containing the data.
#' @param model_obj A `model` object specifying the finite mixture model.
#' @return None. If any check is not passed an error is thrown.
#' @noRd
#' @keywords internal
#' 
#' @seealso 
#' * [fdata][fdata_class] for the definition of the `fdata` class
#' * [model][model_class] for the definition of the `model` class
#' * [mcmc][mcmc_class] for the definition of the `mcmc` class
#' * [mixturemcmc()] for the starting MCMC sampling
".check.mcmc.Mcmcstart" <- function(mcmc.obj) {
  if (class(mcmc.obj) != "mcmc") {
    stop(paste("Wrong argument. 'mcmc' must be an object of class ",
      "'mcmc'.",
      sep = ""
    ))
  }
  .valid.MCMC(mcmc.obj)
}
### Logic
### Logic parameters: Generates starting parameters for @dist in
### 'model.obj'. Returns a 'model' object with starting parameters
### in @par.
#' Sets starting parameters for `mcmcstart`
#' 
#' For internal usage only. This function sets the parameters of a finite 
#' mixture model defined in the slots `@@par` and `@@weight`. 
#' 
#' @param fdata_obj An `fdata` object containing the data.
#' @param model_obj A `model` object specifying the finite mixture model.
#' @param mcmc_obj An `mcmc` object specifying the hyper-parameters for MCMC 
#'   sampling.
#' @return A `model` object with starting parameters. 
#' @noRd
#' @keywords internal
#' 
#' @seealso 
#' * [fdata][fdata_class] for the definition of the `fdata` class
#' * [model][model_class] for the definition of the `model` class
#' * [mcmc][mcmc_class] for the definition of the `mcmc` class
#' * [mixturemcmc()] for the starting MCMC sampling
".parameters.Mcmcstart" <- function(fdata.obj, model.obj, mcmc.obj) {
  K <- model.obj@K
  dist <- model.obj@dist
  ## Check if model object for student-t distributions has
  ## a parameter 'df'.
  if (dist %in% c("student", "studmult")) {
    .mcmcstart.Student.Df(model.obj)
  }
  ## Check if weights have been already initialized
  if (K > 1) {
    if (model.obj@indicmod == "multinomial") {
      model.obj <- .parameters.multinomial.Mcmcstart(model.obj)
    } ## else: Markov model, implemented later.
  }
  if (dist %in% c("poisson", "cond.poisson")) {
    .parameters.poisson.Mcmcstart(fdata.obj, model.obj)
  } else if (dist == "exponential") {
    .parameters.exponential.Mcmcstart(fdata.obj, model.obj, mcmc.obj)
  } else if (dist == "binomial") {
    .parameters.binomial.Mcmcstart(fdata.obj, model.obj)
  } else if (dist == "normal" || dist == "student") {
    .parameters.Norstud.Mcmcstart(fdata.obj, model.obj, mcmc.obj)
  } else if (dist %in% c("normult", "studmult")) {
    .parameters.Norstudmult.Mcmcstart(fdata.obj, model.obj, mcmc.obj)
  }
}

#' Set up exposures in `mcmcstart`
#' 
#' For internal usage only. This function sets up the exposures Poisson mixture 
#' model. If the `fdata` object already contains `exposures` these are checked 
#' for consistency with the number of observations `N`. if exposures cannot be 
#' set an error is thrown.
#' 
#' @param fdata_obj An `fdata` object containing the data.
#' @return A matrix containing th exposures.
#' @noRd
#' @keywords internal
#' 
#' @seealso 
#' * [fdata][fdata_class] for the definition of the `fdata` class
#' * [model][model_class] for the definition of the `model` class
#' * [mcmc][mcmc_class] for the definition of the `mcmc` class
#' * [mixturemcmc()] for the starting MCMC sampling
".mcmcstart.Exp" <- function(data.obj) {
  r <- data.obj@r
  N <- data.obj@N
  has.exp <- (length(data.obj@exp) > 0)
  if (has.exp) {
    if (data.obj@bycolumn) {
      if (nrow(data.obj@exp) != N && nrow(data.obj@exp) != 1) {
        stop(paste(
          "Dimension of slot 'exp' of 'data' object",
          "does not match dimension of slot 'y' of",
          "'data' object."
        ),
        sep = ""
        )
      } else if (nrow(data.obj@exp) == N) {
        exp <- data.obj@exp
      } else {
        exp <- matrix(data.obj@exp[1, 1], nrow = N, ncol = 1)
      }
    } else {
      if (ncol(data.obj@exp) != N && ncol(data.obj@exp) != 1) {
        stop(paste(
          "Dimension of slot 'exp' of 'data' object",
          "does not match dimension of slot 'y' of",
          "'data' object."
        ),
        sep = ""
        )
      } else if (ncol(data.obj@exp) == N) {
        exp <- t(data.obj@exp)
      } else {
        exp <- matrix(data.obj@exp[1, 1], nrow = N, ncol = 1)
      }
    }
  } else {
    exp <- matrix(1, nrow = N, ncol = 1)
  }
  return(exp)
}

#' Set up starting parameters for the weights in `mcmcstart`
#' 
#' For internal usage only. This function sets up the starting weights for a 
#' finite mixture model, by referring to multinomial model for the indicators. 
#' Starting weights are chosen to be equally weighted by the number of 
#' components in the `model` object's slot `@@K`.
#' 
#' @param model_obj A `model` object specifying the finite mixture model.
#' @return A `model` object with starting weights. 
#' @noRd
#' @keywords internal
#' 
#' @seealso 
#' * [fdata][fdata_class] for the definition of the `fdata` class
#' * [model][model_class] for the definition of the `model` class
#' * [mcmc][mcmc_class] for the definition of the `mcmc` class
#' * [mixturemcmc()] for the starting MCMC sampling
".parameters.multinomial.Mcmcstart" <- function(model.obj) {
  K <- model.obj@K
  if (!hasWeight(model.obj)) {
    model.obj@weight <- matrix(1 / K, nrow = 1, ncol = K)
  }
  return(model.obj)
}

#' Set up starting parameters for a Poisson mixture in `mcmcstart`
#' 
#' For internal usage only. This function sets up the starting parameters for 
#' a Poisson mixture model specified by its argument. 
#' 
#' @param fdata_obj An `fdata_obj` storing the data.
#' @param model_obj A `model` object specifying the finite mixture model.
#' @return A `model` object with starting parameters.
#' @noRd
#' @keywords internal
#' 
#' @seealso 
#' * [fdata][fdata_class] for the definition of the `fdata` class
#' * [model][model_class] for the definition of the `model` class
#' * [mcmc][mcmc_class] for the definition of the `mcmc` class
".parameters.poisson.Mcmcstart" <- function(fdata.obj, model.obj) {
  K <- model.obj@K
  datam <- getColY(fdata.obj)
  if (!hasPar(model.obj)) {
    if (K == 1) {
      pm <- max(mean(datam / exp, na.rm = TRUE), 0.1)
      pm <- array(pm, dim = c(1, K))
    } else { ## K > 1
      if (hasExp(fdata.obj)) {
        expos <- getColExp(fdata.obj)
        pm <- (mean(datam / expos, na.rm = TRUE)) * exp(runif(K))
      } else {
        pm <- (mean(datam, na.rm = TRUE)) * exp(runif(K))
      }
      pm <- pmax(pm, 0.1)
    }
    model.obj@par <- list(lambda = pm)
  }
  return(model.obj)
}

#' Set up starting parameters for an exponential mixture in `mcmcstart`
#' 
#' For internal usage only. This function sets up the starting parameters for 
#' an exponential mixture model specified by its argument. 
#' 
#' @param fdata_obj An `fdata_obj` storing the data.
#' @param model_obj A `model` object specifying the finite mixture model.
#' @param mcmc_obj An `mcmc` object containing all hyper-parameters for MCMC
#'   sampling.
#' @return A `model` object with starting parameters.
#' @noRd
#' @keywords internal
#' 
#' @seealso 
#' * [fdata][fdata_class] for the definition of the `fdata` class
#' * [model][model_class] for the definition of the `model` class
#' * [mcmc][mcmc_class] for the definition of the `mcmc` class
".parameters.exponential.Mcmcstart" <- function(fdata.obj, model.obj,
                                                mcmc.obj) {
  if (!hasPar(model.obj)) {
    datam <- getColY(fdata.obj)
    K <- model.obj@K
    if (K == 1) {
      pm <- 1 / mean(datam, na.rm = TRUE)
    } else { ## K > 1
      pm <- exp(runif(K)) / mean(datam, na.rm = TRUE)
    }
    model.obj@par <- list(lambda = pm)
  }
  return(model.obj)
}

#' Set up starting parameters for a Binomial mixture in `mcmcstart`
#' 
#' For internal usage only. This function sets up the starting parameters for 
#' a Binomial mixture model specified by its argument. 
#' 
#' @param fdata_obj An `fdata_obj` storing the data.
#' @param model_obj A `model` object specifying the finite mixture model.
#' @return A `model` object with starting parameters.
#' @noRd
#' @keywords internal
#' 
#' @seealso 
#' * [fdata][fdata_class] for the definition of the `fdata` class
#' * [model][model_class] for the definition of the `model` class
#' * [mcmc][mcmc_class] for the definition of the `mcmc` class
".parameters.binomial.Mcmcstart" <- function(fdata.obj, model.obj) {
  if (!hasPar(model.obj) && hasT(fdata.obj, verbose = TRUE)) {
    datam <- getColY(fdata.obj)
    K <- model.obj@K
    if (K == 1) {
      pm <- mean(datam / fdata.obj@T, na.rm = TRUE)
      pm <- pmin(pmax(pm, 0.1), 0.9)
    } else { ## K > 1
      pm <- mean(datam / fdata.obj@T, na.rm = TRUE) * exp(.2 * runif(K))
      pm <- pmin(pmax(pm, 0.1), 0.9)
    }
    model.obj@par <- list(p = pm)
  }
  return(model.obj)
}

#' Set up starting parameters for a normal or Student-t mixture in `mcmcstart`
#' 
#' For internal usage only. This function sets up the starting parameters for 
#' a normal or Student-t mixture model specified by its argument. 
#' 
#' @param fdata_obj An `fdata_obj` storing the data.
#' @param model_obj A `model` object specifying the finite mixture model.
#' @param mcmc_obj An `mcmc` object containing all hyper-parameters for MCMC
#'   sampling.
#' @return A `model` object with starting parameters.
#' @noRd
#' @keywords internal
#' 
#' @seealso 
#' * [fdata][fdata_class] for the definition of the `fdata` class
#' * [model][model_class] for the definition of the `model` class
#' * [mcmc][mcmc_class] for the definition of the `mcmc` class
".parameters.Norstud.Mcmcstart" <- function(fdata.obj, model.obj,
                                       mcmc.obj) {
  datam <- getColY(fdata.obj)
  K <- model.obj@K
  has.par <- (length(model.obj@par) > 0)
  start.mu <- FALSE
  start.sigma <- FALSE
  if (!has.par) {
    start.mu <- TRUE
    start.sigma <- TRUE
  } else { ## has already parameters
    start.mu <- !"mu" %in% names(model.obj@par)
    start.sigma <- !"sigma" %in% names(model.obj@par)
  }
  if (start.mu) {
    if (K == 1) {
      pm <- mean(datam, na.rm = TRUE)
    } else { ## K > 1
      pm <- mean(datam, na.rm = TRUE) +
        sd(datam, na.rm = TRUE) * runif(K)
      pm <- matrix(pm, nrow = 1, ncol = K)
    }
    if (start.sigma) {
      model.obj@par <- list(mu = pm)
    } else {
      model.obj@par$mu <- pm
    }
  }
  if (start.sigma) {
    pm <- sd(datam, na.rm = TRUE)
    pm <- matrix(pm, nrow = 1, ncol = K)
    model.obj@par$sigma <- pm
  }

  return(model.obj)
}

#' Set up starting parameters for a multivariate mixture in `mcmcstart`
#' 
#' For internal usage only. This function sets up the starting parameters for 
#' a multivariate normal or Student-t mixture model specified by its argument. 
#' 
#' @param fdata_obj An `fdata_obj` storing the data.
#' @param model_obj A `model` object specifying the finite mixture model.
#' @param mcmc_obj An `mcmc` object containing all hyper-parameters for MCMC
#'   sampling.
#' @return A `model` object with starting parameters.
#' @noRd
#' @keywords internal
#' 
#' @seealso 
#' * [fdata][fdata_class] for the definition of the `fdata` class
#' * [model][model_class] for the definition of the `model` class
#' * [mcmc][mcmc_class] for the definition of the `mcmc` class
".parameters.Norstudmult.Mcmcstart" <- function(fdata.obj, model.obj,
                                           mcmc.obj) {
  K <- model.obj@K
  r <- model.obj@r
  has.par <- (length(model.obj@par) > 0)
  datam <- getColY(fdata.obj)
  ## Check if parameters are already provided ##
  start.mu <- FALSE
  start.sigma <- FALSE
  if (!has.par) {
    start.mu <- TRUE
    start.sigma <- TRUE
  } else {
    has.mu <- "mu" %in% names(model.obj@par)
    has.sigma <- "sigma" %in% names(model.obj@par)
    if (!has.mu) {
      start.mu <- TRUE
    }
    if (!has.sigma) {
      start.sigma <- TRUE
    }
  }
  cov.m <- cov(datam)
  if (start.mu) {
    if (K == 1) {
      pm.mu <- apply(datam, 2, mean, na.rm = TRUE)
      pm.mu <- array(pm.mu, dim = c(1, K))
    } else { ## K > 1
      mean <- apply(datam, 2, mean, na.rm = TRUE)
      pm.mu <- matrix(0, nrow = r, ncol = K)
      for (i in 1:K) {
        pm.mu[, i] <- matrix(mean) + t(chol(cov.m)) %*% matrix(runif(K))
      }
    }
    if (!has.par) {
      model.obj@par <- list(mu = pm.mu)
    } else {
      model.obj@par$mu <- pm.mu
    }
  }
  if (start.sigma) {
    model.obj@par$sigma <- array(cov.m, dim = c(r, r, K))
  }
  return(model.obj)
}

#' Set up starting degrees of freedom for a Student-t mixture in `mcmcstart`
#' 
#' For internal usage only. This function sets up the starting degrees of 
#' freedom for a Student-t mixture model specified by its argument. 
#' 
#' @param model_obj A `model` object specifying the finite mixture model.
#' @return A `model` object with starting parameters.
#' @noRd
#' @keywords internal
#' 
#' @seealso 
#' * [fdata][fdata_class] for the definition of the `fdata` class
#' * [model][model_class] for the definition of the `model` class
#' * [mcmc][mcmc_class] for the definition of the `mcmc` class
".mcmcstart.Student.Df" <- function(model.obj) {
  K <- model.obj@K
  has.par <- (length(model.obj@par) > 0)
  if (has.par) {
    has.df <- "df" %in% names(model.obj@par)
    if (!has.df) {
      model.obj@pari$df <- array(10, dim = c(1, K))
      validObject(model.obj)
    }
  } else {
    model@par <- list(df = array(10, dim = c(1, K)))
  }
  return(model.obj)
}

### Logic indicators: Returns an 'fdata' object with generated
### indicators.
#' Set up starting indicators for a finite mixture in `mcmcstart`
#' 
#' For internal usage only. This function sets up the starting indicators for 
#' a finite mixture model specified by its argument. 
#' 
#' @param fdata_obj An `fdata_obj` storing the data.
#' @param model_obj A `model` object specifying the finite mixture model.
#' @return A `model` object with starting parameters.
#' @noRd
#' @keywords internal
#' 
#' @seealso 
#' * [fdata][fdata_class] for the definition of the `fdata` class
#' * [model][model_class] for the definition of the `model` class
#' * [mcmc][mcmc_class] for the definition of the `mcmc` class
".indicators.Mcmcstart" <- function(fdata.obj, model.obj) {
  dist <- model.obj@dist
  if (dist %in% c("poisson", "cond.poisson", "exponential")) {
    .indicators.poisson.Mcmcstart(fdata.obj, model.obj)
  } else if (dist == "binomial") {
    .indicators.binomial.Mcmcstart(fdata.obj, model.obj)
  } else if (dist %in% c(
    "normal", "normult",
    "student", "studmult"
  )) {
    .mcmcstart.Ind.Norstud(fdata.obj, model.obj)
  }
}

### Logic indicators for Poisson: If it is started by sampling
### the parameters a simple kmeans-clustering is performed
### to find initial indicators. If indicators are already
### in slot @S of the 'fdata' object, the 'fdata' object is
### immediately returned.
#' Set up starting indicators for a Poisson mixture in `mcmcstart`
#' 
#' For internal usage only. This function sets up the starting indicators for 
#' a Poisson or exponential mixture model specified by its argument. 
#' 
#' @param fdata_obj An `fdata_obj` storing the data.
#' @param model_obj A `model` object specifying the finite mixture model.
#' @return A `model` object with starting parameters.
#' @noRd
#' @keywords internal
#' 
#' @seealso 
#' * [fdata][fdata_class] for the definition of the `fdata` class
#' * [model][model_class] for the definition of the `model` class
#' * [mcmc][mcmc_class] for the definition of the `mcmc` class
".indicators.poisson.Mcmcstart" <- function(fdata.obj, model.obj) {
  K <- model.obj@K
  if (!hasS(fdata.obj)) {
    datam <- getColY(fdata.obj)
    S <- matrix(kmeans(datam^.5,
      centers = K,
      nstart = K
    )$cluster)
    if (fdata.obj@bycolumn) {
      fdata.obj@S <- S
    } else {
      fdata.obj@S <- t(S)
    }
  }
  return(fdata.obj)
}

#' Set up starting indicators for a Binomial mixture in `mcmcstart`
#' 
#' For internal usage only. This function sets up the starting indicators for 
#' a Binomial mixture model specified by its argument. 
#' 
#' @param fdata_obj An `fdata_obj` storing the data.
#' @param model_obj A `model` object specifying the finite mixture model.
#' @return A `model` object with starting parameters.
#' @noRd
#' @keywords internal
#' 
#' @seealso 
#' * [fdata][fdata_class] for the definition of the `fdata` class
#' * [model][model_class] for the definition of the `model` class
#' * [mcmc][mcmc_class] for the definition of the `mcmc` class
".indicators.binomial.Mcmcstart" <- function(fdata.obj, model.obj) {
  if (!hasS(fdata.obj)) {
    K <- model.obj@K
    datam <- getColY(fdata.obj)
    if ((max(datam) - min(datam)) > 2 * K) {
      ## use k-means to determine a starting classification
      if (fdata.obj@bycolumn) {
        fdata.obj@S <- as.matrix(kmeans(datam^.5,
          centers = K,
          nstart = K
        )$cluster)
      } else {
        fdata.obj@S <- t(as.matrix(kmeans(datam^.5,
          centers = K,
          nstart = K
        )$cluster))
      }
    } else {
      ## random classification
      N <- fdata.obj@N
      if (fdata.obj@bycolumn) {
        fdata.obj@S <- as.matrix(sample(c(1:K), N,
          replace = TRUE
        ))
      } else {
        fdata.obj@S <- t(as.matrix(sample(c(1:K), N,
          replace = TRUE
        )))
      }
    }
  }
  return(fdata.obj)
}

#' Set up starting indicators for a normal or Student-t mixture in `mcmcstart`
#' 
#' For internal usage only. This function sets up the starting indicators for 
#' a normal or Student-t mixture model specified by its argument. 
#' 
#' @param fdata_obj An `fdata_obj` storing the data.
#' @param model_obj A `model` object specifying the finite mixture model.
#' @return A `model` object with starting parameters.
#' @noRd
#' @keywords internal
#' 
#' @seealso 
#' * [fdata][fdata_class] for the definition of the `fdata` class
#' * [model][model_class] for the definition of the `model` class
#' * [mcmc][mcmc_class] for the definition of the `mcmc` class
".mcmcstart.Ind.Norstud" <- function(data.obj, model.obj) {
  K <- model.obj@K
  # Checks, if slot 'S' ist in 'data.obj'. If not, throws an error
  .valid.S.Fdata(data.obj)
  # Because of the line above, we can set this.``
  has.S <- TRUE
  if (has.S) {
    return(data.obj)
  } else {
    if (data.obj@bycolumn) {
      .valid.y.Fdata(data.obj)
      datam <- data.obj@y
      data.obj@S <- as.matrix(kmeans(datam^.5,
        centers = K,
        nstart = K
      )$cluster)
    } else {
      data.obj@S <- t(as.matrix(kmeans(datam^.5,
        centers = K,
        nstart = K
      )$cluster))
    }
    return(data.obj)
  }
}
