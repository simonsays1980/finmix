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
# USED IN modelmoments 

#' @noRd
".mixturemoments.normal" <- function(model, J, meanm) {
  zm <- array(0, dim = c(J, 1))
  zm[seq(2, J, by = 2)] <- exp(cumsum(log(seq(1, (J - 1), by = 2))))
  moments <- array(0, dim = c(J, 1))

  for (m in 2:J) { ## first higher moment is always zero
    diff <- model@par$mu - meanm
    moments[m] <- sum(model@weight * diff^m)
    for (n in 1:m) {
      cm <- diff^(m - n) * model@par$sigma^(n / 2) * zm[n]
      moments[m] <- moments[m] + choose(m, n) * sum(model@weight * cm)
    }
  } 

  return(moments)
}

#' @noRd
".mixturemoments.student" <- function(model, J, meanm) {
  moments <- array(0, dim = c(J, 1))
  sigma <- model@par$sigma
  mu <- model@par$mu
  degrees <- model@par$df
  weight <- model@weight
  diff <- mu - meanm
  raw.moments <- array(0, dim = c(1, model@K))
  for (j in seq(1, J)) {
    moments[j] <- sum(diff^j * weight)
    for (n in seq(1, j)) {
      raw.moments <- .raw.moments.student(n, sigma, degrees)
      moments[j] <- moments[j] + sum(choose(j, n) *
        diff^(j - n) * raw.moments * weight)
    }
  }
  return(moments)
}

#' @noRd
".mixturemoments.exponential" <- function(model, J, meanm) {
  moments <- array(0, dim = c(J, 1))
  lambda <- model@par$lambda
  weight <- model@weight
  diff <- 1 / lambda - meanm
  for (j in seq(1, J)) {
    moments[j] <- sum(diff^j * weight) ## case E[(X-mu)^0]
    for (n in seq(1, j)) {
      raw.moments <- .raw.moments.exponential(n, lambda)
      moments[j] <- moments[j] + sum(choose(j, n) *
        diff^(j - n) * raw.moments * weight)
    }
  }
  return(moments)
}

#' @noRd
".raw.moments.student" <- function(n, sigma, degrees) {
  value <- array(0, dim = c(1, length(degrees)))
  if (n > 0 && n %% 2 == 0) {
    for (i in seq(1, n / 2)) {
      value <- (2 * i - 1) / (degrees - 2 * i) * degrees^(n / 2) * sqrt(sigma)^n
    }
    value[degrees <= n] <- Inf
  } else {
    value[degrees <= n] <- NaN
  }
  return(value)
}

#' @noRd
".raw.moments.exponential" <- function(n, lambda) {
  values <- rep(0, length(lambda))
  for (i in seq(0, n)) {
    values <- values + factorial(n) / lambda^n * (-1)^i / factorial(i)
  }
  return(values)
}