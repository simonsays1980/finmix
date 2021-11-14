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

#' Computes the log-likelihood for normal mixture components
#' 
#' @description
#' For internal usage only. This function calculates the likelihood of each 
#' mixture component for the data. In addition the maximum likelihood along all 
#' mixture components and the log-likelihood is calculated. 
#' 
#' @param y A matrix containing the data. Of dimension `Nx1` for 
#'   univariate models and `Nxr` for multivariate models.
#' @param mu A vector containing the means of the normal mixture components.
#' @param sigma A vector containing the standard deviations of the normal 
#'   mixture components.
#' @return A list containing the likelihood, the maximum likelihood and the 
#'   log-likelihood.
#' @noRd
#' 
#' @seealso 
#' * [fdata-class] for the `fdata` class definition
#' * [model][model_class] for the `model` class definition
#' * [dataclass()] for the constructor of the `dataclass` class
".likelihood.normal" <- function(y, mu, sigma) {
  N <- nrow(y)
  K <- ncol(mu)
  y <- matrix(y, nrow = N, ncol = K)

  err <- t(apply(y, 1, "-", mu))
  err <- t(apply(err^2, 1, "/", sigma))

  loglik <- -.5 * (log(2 * pi) + t(apply(err, 1, "+", log(sigma))))

  if (K == 1) {
    max.lik <- loglik
  } else {
    max.lik <- apply(loglik, 1, max, na.rm = TRUE)
  }
  l.h <- exp(loglik - max.lik)

  result <- list(lh = l.h, maxl = matrix(max.lik), llh = loglik)
  return(result)
}

#' Computes the log-likelihood for student mixture components
#' 
#' @description
#' For internal usage only. This function calculates the likelihood of each 
#' mixture component for the data. In addition the maximum likelihood along all 
#' mixture components and the log-likelihood is calculated. 
#' 
#' @param y A matrix containing the data. Of dimension `Nx1` for 
#'   univariate models and `Nxr` for multivariate models.
#' @param mu A vector containing the means of the student mixture components.
#' @param sigma A vector containing the standard deviations of the student 
#'   mixture components.
#' @param df A vector containing the degrees of freedom of the student mixture 
#'   components.
#' @return A list containing the likelihood, the maximum likelihood and the 
#'   log-likelihood.
#' @noRd
#' 
#' @seealso 
#' * [fdata-class] for the `fdata` class definition
#' * [model][model_class] for the `model` class definition
#' * [dataclass()] for the constructor of the `dataclass` class
".likelihood.student" <- function(y, mu, sigma, df) {
  N <- nrow(y)
  K <- ncol(mu)
  y <- matrix(y, nrow = N, ncol = K)
  mu <- matrix(mu, nrow = N, ncol = K, byrow = TRUE)
  sigma <- matrix(sigma, nrow = N, ncol = K, byrow = TRUE)
  df <- matrix(df, nrow = N, ncol = K, byrow = TRUE)

  err <- (y - mu)^2 / sigma


  loglik <- lgamma((df + 1) / 2) - lgamma(df / 2) - .5 * (log(df * pi) + log(sigma)) - (df + 1) / 2 * log(1 + err / df)

  if (K == 1) {
    max.lik <- loglik
  } else {
    max.lik <- apply(loglik, 1, max, na.rm = TRUE)
  }
  l.h <- exp(loglik - max.lik)

  result <- list(lh = l.h, maxl = matrix(max.lik), llh = loglik)
  return(result)
}

#' Computes the log-likelihood for exponential mixture components
#' 
#' @description
#' For internal usage only. This function calculates the likelihood of each 
#' mixture component for the data. In addition the maximum likelihood along all 
#' mixture components and the log-likelihood is calculated. 
#' 
#' @param y A matrix containing the data. Of dimension `Nx1` for 
#'   univariate models and `Nxr` for multivariate models.
#' @param lambda A vector containing the rates of the exponential mixture 
#'   components.
#' @return A list containing the likelihood, the maximum likelihood and the 
#'   log-likelihood.
#' @noRd
#' 
#' @seealso 
#' * [fdata-class] for the `fdata` class definition
#' * [model][model_class] for the `model` class definition
#' * [dataclass()] for the constructor of the `dataclass` class
".likelihood.exponential" <- function(y, lambda) {
  N <- nrow(y)
  K <- ncol(lambda)
  y <- matrix(y, nrow = N, ncol = K)

  lambda <- matrix(lambda, nrow = N, ncol = K, byrow = TRUE)

  lambda <- apply(lambda, c(1, 2), max, 0.0001)
  loglik <- log(lambda) - y * lambda
  max.lik <- apply(loglik, 1, max)

  l.h <- exp(loglik - max.lik)

  result <- list(lh = l.h, maxl = matrix(max.lik), llh = loglik)
  return(result)
}

#' Computes the log-likelihood for Poisson mixture components
#' 
#' @description
#' For internal usage only. This function calculates the likelihood of each 
#' mixture component for the data. In addition the maximum likelihood along all 
#' mixture components and the log-likelihood is calculated. 
#' 
#' @param y A matrix containing the data. Of dimension `Nx1` for 
#'   univariate models and `Nxr` for multivariate models.
#' @param lambda A vector containing the rates of the Poisson mixture 
#'   components.
#' @return A list containing the likelihood, the maximum likelihood and the 
#'   log-likelihood.
#' @noRd
#' 
#' @seealso 
#' * [fdata-class] for the `fdata` class definition
#' * [model][model_class] for the `model` class definition
#' * [dataclass()] for the constructor of the `dataclass` class
".likelihood.poisson" <- function(y, lambda) {
  N <- nrow(y)
  K <- ncol(lambda)
  nst <- nrow(lambda)

  y <- matrix(y, nrow = N, ncol = K)

  if (nst == 1) {
    lambda <- matrix(lambda, nrow = N, ncol = K, byrow = TRUE)
  }
  lambda <- apply(lambda, c(1, 2), max, 10e-5, na.rm = TRUE)

  loglik <- y * log(lambda) - lambda - lgamma(y + 1)

  max.lik <- apply(loglik, 1, max, na.rm = TRUE)

  l.h <- exp(apply(loglik, 2, "-", max.lik))

  result <- list(lh = l.h, maxl = matrix(max.lik), llh = loglik)
  return(result)
}

#' Computes the log-likelihood for Binomial mixture components
#' 
#' @description
#' For internal usage only. This function calculates the likelihood of each 
#' mixture component for the data. In addition the maximum likelihood along all 
#' mixture components and the log-likelihood is calculated. 
#' 
#' @param y A matrix containing the data. Of dimension `Nx1` for 
#'   univariate models and `Nxr` for multivariate models.
#' @param lambda A vector containing the probabilities of the Binomial mixture 
#'   components.
#' @return A list containing the likelihood, the maximum likelihood and the 
#'   log-likelihood.
#' @noRd
#' 
#' @seealso 
#' * [fdata-class] for the `fdata` class definition
#' * [model][model_class] for the `model` class definition
#' * [dataclass()] for the constructor of the `dataclass` class
".likelihood.binomial" <- function(y, T, p) {
  N <- nrow(y)
  K <- length(p)
  nst <- nrow(T)

  y <- matrix(y, nrow = N, ncol = K)
  T <- matrix(T, nrow = N, ncol = K, byrow = TRUE)

  loglik <- lgamma(T + 1) - lgamma(T - y + 1) - lgamma(y + 1)
  loglik <- loglik + t((apply(y, 1, "*", log(p)))) + t((apply(T - y, 1, "*", log(1 - p))))

  max.lik <- apply(loglik, 1, max, na.rm = TRUE)

  l.h <- exp(apply(loglik, 2, "-", max.lik))

  result <- list(lh = l.h, maxl = matrix(max.lik), llh = loglik)
  return(result)
}

#' Computes the log-likelihood for multivariate normal mixture components
#' 
#' @description
#' For internal usage only. This function calculates the likelihood of each 
#' mixture component for the data. In addition the maximum likelihood along all 
#' mixture components and the log-likelihood is calculated. 
#' 
#' @param y A matrix containing the data. Of dimension `Nx1` for 
#'   univariate models and `Nxr` for multivariate models.
#' @param mu A matrix containing the means for each mixture component. Of 
#'   dimension `rxK`
#' @param sigmainv An array containing the inverse variance-covariance matrices 
#'   for each mixture component. Of dimension `rxrxK`.
#' @param logdet A vector containing the logarithmized determinants of the 
#'   variance-covariance matrices for each component. Of dimension `1xK`.
#' @return A list containing the likelihood, the maximum likelihood and the 
#'   log-likelihood.
#' @noRd
#' 
#' @seealso 
#' * [fdata-class] for the `fdata` class definition
#' * [model][model_class] for the `model` class definition
#' * [dataclass()] for the constructor of the `dataclass` class
".likelihood.normult" <- function(y, mu, sigmainv, logdet) {
  N <- nrow(y)
  r <- ncol(y)
  K <- dim(sigmainv)[3]
  loglik <- matrix(0, nrow = N, ncol = K)
  loglik1 <- -.5 * r * log(2 * pi)

  for (k in 1:K) {
    eps <- t(apply(y, 1, "-", t(matrix(mu[, k]))))
    loglik[, k] <- loglik1 + .5 * logdet[k] - .5 * apply(eps %*% sigmainv[, , k] * eps, 1, sum)
  }

  maxlik <- t(apply(loglik, 1, max))
  l.h <- exp(apply(loglik, 2, "-", maxlik))

  results <- list(lh = l.h, maxl = matrix(maxlik), llh = loglik)
  return(results)
}

#' Computes the log-likelihood for multivariate Student-t mixture components
#' 
#' @description
#' For internal usage only. This function calculates the likelihood of each 
#' mixture component for the data. In addition the maximum likelihood along all 
#' mixture components and the log-likelihood is calculated. 
#' 
#' @param y A matrix containing the data. Of dimension `Nx1` for 
#'   univariate models and `Nxr` for multivariate models.
#' @param mu A matrix containing the means for each mixture component. Of 
#'   dimension `rxK`
#' @param sigmainv An array containing the inverse variance-covariance matrices 
#'   for each mixture component. Of dimension `rxrxK`.
#' @param logdet A vector containing the logarithmized determinants of the 
#'   variance-covariance matrices for each component. Of dimension `1xK`.
#' @return A list containing the likelihood, the maximum likelihood and the 
#'   log-likelihood.
#' @noRd
#' 
#' @seealso 
#' * [fdata-class] for the `fdata` class definition
#' * [model][model_class] for the `model` class definition
#' * [dataclass()] for the constructor of the `dataclass` class
".likelihood.studmult" <- function(y, mu, sigmainv, logdet, df) {
  N <- nrow(y)
  K <- ncol(mu)
  r <- ncol(y)

  loglik <- matrix(0, nrow = N, ncol = K)

  for (k in 1:K) {
    mum <- matrix(mu[, k], nrow = N, ncol = r, byrow = TRUE)
    err <- y - mum
    err <- err %*% sigmainv[, , k] * err
    err <- apply(err, 1, sum)
    loglik[, k] <- lgamma((df[k] + r) / 2) - lgamma(df[k] / 2) + .5 * logdet[k] - .5 * r * log(df[k] * pi)
    loglik[, k] <- loglik[, k] - ((df[k] + r) / 2) * log(1 + err / df[k])
  }

  maxlik <- t(apply(loglik, 1, max))
  l.h <- exp(apply(loglik, 2, "-", maxlik))

  result <- list(lh = l.h, maxl = matrix(maxlik), llh = loglik)
  return(result)
}
