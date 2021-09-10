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

.mcmcoutputfixpost <- setClass("mcmcoutputfixpost",
  representation(post = "list"),
  contains = c("mcmcoutputfix"),
  validity = function(object) {
    ## else: OK
    TRUE
  },
  prototype(post = list())
)

setMethod(
  "show", "mcmcoutputfixpost",
  function(object) {
    cat("Object 'mcmcoutputfixpost\n")
    cat("     class       :", class(object), "\n")
    cat("     M           :", object@M, "\n")
    cat("     burnin      :", object@burnin, "\n")
    cat("     ranperm     :", object@ranperm, "\n")
    cat(
      "     par         : List of",
      length(object@par), "\n"
    )
    cat(
      "     log         : List of",
      length(object@log), "\n"
    )
    cat(
      "     post        : List of",
      length(object@post), "\n"
    )
    cat(
      "     model       : Object of class",
      class(object@model), "\n"
    )
    cat(
      "     prior       : Object of class",
      class(object@prior), "\n"
    )
  }
)

setMethod(
  "plotTraces", signature(
    x = "mcmcoutputfixpost",
    dev = "ANY",
    lik = "ANY",
    col = "ANY"
  ),
  function(x, dev = TRUE, lik = 1, col = FALSE, ...) {
    ## Call 'plot()' from 'mcmcoutputfix
    callNextMethod(x, dev, lik, col, ...)
  }
)

setMethod(
  "plotHist", signature(
    x = "mcmcoutputfixpost",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    ## Call 'plotHist()' from 'mcmcoutputfix'
    callNextMethod(x, dev, ...)
  }
)

setMethod(
  "plotDens", signature(
    x = "mcmcoutputfixpost",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    ## Call 'plotDens()' from 'mcmcoutputfix'
    callNextMethod(x, dev, ...)
  }
)

setMethod(
  "plotPointProc", signature(
    x = "mcmcoutputfixpost",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    ## Call 'plotPointProc()' from 'mcmcoutputfix'
    callNextMethod(x, dev, ...)
  }
)

setMethod(
  "plotSampRep", signature(
    x = "mcmcoutputfixpost",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    ## Call 'plotSampRep()' from 'mcmcoutputfix'
    callNextMethod(x, dev, ...)
  }
)

setMethod(
  "plotPostDens", signature(
    x = "mcmcoutputfixpost",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    ## Call 'plotPostDens()' from 'mcmcoutputfix'
    callNextMethod(x, dev, ...)
  }
)

setMethod(
  "subseq", signature(
    object = "mcmcoutputfixpost",
    index = "array"
  ),
  function(object, index) {
    ## Call 'subseq()' from 'mcmcoutputfix'
    callNextMethod(object, index)
    dist <- object@model@dist
    ## post ##
    if (dist == "poisson") {
      .subseq.Poisson.Post(object, index)
    } else if (dist == "binomial") {
      .subseq.Binomial.Mcmcoutputfixpost(object, index)
    } else if (dist %in% c("normal", "student")) {
      .subseq.Norstud.Mcmcoutputfixpost(object, index)
    } else if (dist %in% c("normult", "studmult")) {
      .subseq.Normultstud.Mcmcoutputfixpost(object, index)
    }
  }
)

setMethod(
  "swapElements", signature(
    object = "mcmcoutputfixpost",
    index = "array"
  ),
  function(object, index) {
    if (object@model@K == 1) {
      return(object)
    } else {
      ## Call method 'swapiElements()' from 'mcmcoutputfix'
      object <- callNextMethod()
      dist <- object@model@dist
      if (dist == "poisson") {
        .swapElements.Poisson(object, index)
      } else if (dist == "binomial") {
        .swapElements.Binomial.Mcmcoutputfixpost(object, index)
      } else if (dist %in% c("normal", "student")) {
        .swapElements.Norstud.Mcmcoutputfixpost(object, index)
      } else if (dist %in% c("normult", "studmult")) {
        .swapElements.Normultstud.Mcmcoutputfixpost(object, index)
      }
    }
  }
)

setMethod(
  "getPost", "mcmcoutputfixpost",
  function(object) {
    return(object@post)
  }
)

## No setters as users are not intended to manipulate ##
## this object. ##

".subseq.Poisson.Post" <- function(obj, index) {
  if (obj@model@K == 1) {
    obj@post$par$a <- array(obj@post$par$a[index],
      dim = c(obj@M, 1)
    )
    obj@post$par$b <- array(obj@post$par$b[index],
      dim = c(obj@M, 1)
    )
  } else {
    obj@post$par$a <- obj@post$par$a[index, ]
    obj@post$par$b <- obj@post$par$b[index, ]
  }
  return(obj)
}

".swapElements.Poisson.Post" <- function(obj, index) {
  ## Rcpp::export 'swap_cc'
  obj@post$par$a <- swap_cc(obj@post$par$a, index)
  obj@post$par$b <- swap_cc(obj@post$par$b, index)
  return(obj)
}

".subseq.Binomial.Mcmcoutputfixpost" <- function(obj, index) {
  if (obj@model@K == 1) {
    obj@post$par$a <- array(obj@post$par$a[index],
      dim = c(obj@M, 1)
    )
    obj@post$par$b <- array(obj@post$par$b[index],
      dim = c(obj@M, 1)
    )
  } else {
    obj@post$par$a <- obj@post$par$a[index, ]
    obj@post$par$b <- obj@post$par$b[index, ]
  }
  return(obj)
}

".swapElements.Binomial.Mcmcoutputfixpost" <- function(obj, index) {
  ## Rcpp::export 'swap_cc'
  obj@post$par$a <- swap_cc(obj@post$par$a, index)
  obj@post$par$b <- swap_cc(obj@post$par$b, index)
  return(obj)
}

".subseq.Norstud.Mcmcoutputfixpost" <- function(obj, index) {
  if (obj@model@K == 1) {
    obj@post$par$mu$b <- array(obj@post$par$mu$b[index],
      dim = c(obj@M, 1)
    )
    obj@post$par$mu$B <- array(obj@post$par$mu$B[index],
      dim = c(obj@M, 1)
    )
    obj@post$par$sigma$c <- array(obj@post$par$sigma$c[index],
      dim = c(obj@M, 1)
    )
    obj@post$par$sigma$C <- array(obj@post$par$sigma$C[index],
      dim = c(obj@M, 1)
    )
  } else {
    obj@post$par$mu$b <- obj@post$par$mu$b[index, ]
    obj@post$par$mu$b <- obj@post$par$mu$B[index, ]
    obj@post$par$sigma$c <- obj@post$par$sigma$c[index, ]
    obj@post$par$sigma$C <- obj@post$par$sigma$C[index, ]
  }
  return(obj)
}

".swapElements.Norstud.Mcmcoutputfixpost" <- function(obj, index) {
  ## Rcpp::export 'swap_cc'
  obj@post$par$mu$b <- swap_cc(obj@post$par$mu$b, index)
  obj@post$par$mu$B <- swap_cc(obj@post$par$mu$B, index)
  obj@post$par$sigma$c <- swap_cc(obj@post$par$sigma$c, index)
  obj@post$par$sigma$C <- swap_cc(obj@post$par$sigma$C, index)
  return(obj)
}

".subseq.Normultstud.Mcmcoutputfixpost" <- function(obj, index) {
  if (obj@model@K == 1) {
    obj@post$par$mu$b <- obj@post$par$mu$b[index, ]
    obj@post$par$mu$B <- obj@post$par$mu$B[index, ]
    obj@post$par$sigma$c <- obj@post$par$sigma$c[index, ]
    obj@post$par$sigma$C <- obj@post$par$sigma$C[index, ]
  } else {
    obj@post$par$mu$b <- obj@post$par$mu$b[index, , ]
    obj@post$par$mu$B <- obj@post$par$mu$B[index, , ]
    obj@post$par$sigma$c <- obj@post$par$sigma$c[index, , ]
    obj@post$par$sigma$C <- obj@post$par$sigma$C[index, , ]
  }
  return(obj)
}

".swapElements.Normultstud.Mcmcoutputfixpost" <- function(obj, index) {
  ## Rcpp::export 'swap_3d_cc'
  obj@post$par$mu$b <- swap_cc(obj@post$par$mu$b, index)
  obj@post$par$mu$B <- swap_3d_cc(obj@post$par$mu$B, index)
  obj@post$par$sigma$c <- swap_cc(obj@post$par$sigma$c, index)
  obj@post$par$sigma$C <- swap_3d_cc(obj@post$par$sigma$C, index)
  return(obj)
}
