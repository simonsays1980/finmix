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

#' Density function of a Student-t distribution
#' 
#' @description 
#' Unused at this moment.
#' 
#' @param x A vector of valued for which the density should be calculated.
#' @param mu A vector containing the mean of the distribution.
#' @param sigma A vector containing the standard deviation of the distribution.
#' @param df A vector containing the degrees of freedom of the distribution.
#' @return The density of the Student-t distribution for the values of `x`. 
#' @keywords internal
"dstud" <- function(x, mu, sigma, df) {
  fun <- gamma((df + 1) / 2) / (gamma(df / 2) * sqrt(df * pi * sigma)) * 
    (1 + (x - mu)^2 / (df * sigma))^(-(df + 1) / 2)

  return(fun)
}
