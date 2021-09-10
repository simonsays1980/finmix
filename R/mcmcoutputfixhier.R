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

.mcmcoutputfixhier <- setClass("mcmcoutputfixhier",
  representation(hyper = "list"),
  contains = c("mcmcoutputfix"),
  validity = function(object) {
    ## else: OK
    TRUE
  },
  prototype(hyper = list())
)

setMethod(
  "show", "mcmcoutputfixhier",
  function(object) {
    cat("Object 'mcmcoutput'\n")
    cat(
      "     class       :", class(object),
      "\n"
    )
    cat("     M           :", object@M, "\n")
    cat("     burnin      :", object@burnin, "\n")
    cat(
      "     ranperm     :", object@ranperm,
      "\n"
    )
    cat(
      "     par         : List of",
      length(object@par), "\n"
    )
    cat(
      "     log         : List of",
      length(object@log), "\n"
    )
    cat(
      "     hyper       : List of",
      length(object@hyper), "\n"
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
    x = "mcmcoutputfixhier",
    dev = "ANY",
    lik = "ANY",
    col = "ANY"
  ),
  function(x, dev = TRUE, lik = 1, col = FALSE, ...) {
    dist <- x@model@dist
    if (lik %in% c(0, 1)) {
      if (dist == "poisson") {
        .traces.Poisson.Hier(x, dev)
      } else if (dist == "binomial") {
        .traces.Binomial(x, dev)
      } else if (dist == "exponential") {
        callNextMethod(x, dev)
      } else if (dist == "normal") {
        .traces.Normal.Hier(x, dev)
      } else if (dist == "student") {
        .traces.Student.Hier(x, dev)
      } else if (dist == "normult") {
        .traces.Normult.Hier(x, dev, col)
      } else if (dist == "studmult") {
        .traces.Studmult.Hier(x, dev, col)
      }
    }
    if (lik %in% c(1, 2)) {
      ## log ##
      .traces.Log(x, dev, col)
    }
  }
)

setMethod(
  "plotHist", signature(
    x = "mcmcoutputfixhier",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    dist <- x@model@dist
    if (dist == "poisson") {
      .hist.Poisson.Hier(x, dev)
    } else if (dist == "binomial") {
      .hist.Binomial(x, dev)
    } else if (dist == "exponential") {
      .hist.Exponential(x, dev)
    } else if (dist == "normal") {
      .hist.Normal.Hier(x, dev)
    } else if (dist == "student") {
      .hist.Student.Hier(x, dev)
    } else if (dist == "normult") {
      .hist.Normult.Hier(x, dev)
    } else if (dist == "studmult") {
      .hist.Studmult.Hier(x, dev)
    }
  }
)

setMethod(
  "plotDens", signature(
    x = "mcmcoutputfixhier",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    dist <- x@model@dist
    if (dist == "poisson") {
      .dens.Poisson.Hier(x, dev)
    } else if (dist == "binomial") {
      .dens.Binomial(x, dev)
    } else if (dist == "exponential") {
      .dens.Exponential(x, dev)
    } else if (dist == "normal") {
      .dens.Normal.Hier(x, dev)
    } else if (dist == "student") {
      .dens.Student.Hier(x, dev)
    } else if (dist == "normult") {
      .dens.Normult.Hier(x, dev)
    } else if (dist == "studmult") {
      .dens.Studmult.Hier(x, dev)
    }
  }
)

setMethod(
  "plotPointProc", signature(
    x = "mcmcoutputfixhier",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    ## Call 'plotPointProc()' from 'mcmcoutputfix'
    callNextMethod(x, dev, ...)
  }
)

setMethod(
  "plotSampRep", signature(
    x = "mcmcoutputfixhier",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    ## Call 'plotSampRep()' from 'mcmcoutputfix'
    callNextMethod(x, dev, ...)
  }
)

setMethod(
  "plotPostDens", signature(
    x = "mcmcoutputfixhier",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    ## Call 'plotPostDens()' from 'mcmcoutputfix'
    callNextMethod(x, dev, ...)
  }
)

setMethod(
  "subseq", signature(
    object = "mcmcoutputfixhier",
    index = "array"
  ),
  function(object, index) {
    ## Call 'subseq()' from 'mcmcoutputfix'
    callNextMethod(object, index)
    dist <- object@model@dist
    ## hyper ##
    if (dist == "poisson") {
      .subseq.Poisson.Hier(object, index)
    } else if (dist == "normal" || dist == "student") {
      .subseq.Normal.Hier(object, index)
    } else if (dist %in% c("normal", "student")) {
      .subseq.Norstud.Hier.(object, index)
    } else if (dist %in% c("normult", "studmult")) {
      .subseq.Normultstud.Hier(object, index)
    }
  }
)

setMethod(
  "swapElements", signature(
    object = "mcmcoutputfixhier",
    index = "array"
  ),
  function(object, index) {
    ## Check arguments, TODO: .validObject ##
    .swapElements.valid.Arg(object, index)
    if (object@model@K == 1) {
      return(object)
    } else {
      ## Call method 'swap()' from 'mcmcoutputfix'
      callNextMethod(object, index)
    }
  }
)

setMethod(
  "getHyper", "mcmcoutputfixhier",
  function(object) {
    return(object@hyper)
  }
)

## No setters for this object as it is not intended ##
## that users manipulate this object. 		    	##

### Private functions
### These functions are not exported.

### Plot

### Plot Traces
### Plot traces Poisson: Plots traces for each component
### parameter of a Poisson mixture and the hyper parameter 'b'.
".traces.Poisson.Hier" <- function(x, dev) {
  K <- x@model@K
  trace.n <- K + 1
  if (.check.grDevice() && y) {
    dev.new(title = "Traceplots")
  }
  par(
    mfrow = c(trace.n, 1), mar = c(1, 0, 0, 0),
    oma = c(4, 5, 4, 4)
  )
  lambda <- x@par$lambda
  for (k in 1:K) {
    plot(lambda[, k],
      type = "l", axes = F,
      col = "gray20", xlab = "", ylab = ""
    )
    axis(2, las = 2, cex.axis = 0.7)
    mtext(
      side = 2, las = 2, bquote(lambda[k = .(k)]),
      cex = 0.6, line = 3
    )
  }
  b <- x@hyper$b
  plot(b,
    type = "l", axes = F,
    col = "gray68", xlab = "", ylab = ""
  )
  axis(2, las = 2, cex.axis = 0.7)
  mtext(side = 2, las = 2, "b", cex = 0.6, line = 3)
  axis(1)
  mtext(side = 1, "Iterations", cex = 0.7, line = 3)
}

".traces.Normal.Hier" <- function(x, dev) {
  K <- x@model@K
  trace.n <- 2 * K + 1
  if (.check.grDevice() && dev) {
    dev.new(title = "Traceplots")
  }
  par(
    mfrow = c(trace.n, 1), mar = c(1, 0, 0, 0),
    oma = c(4, 5, 4, 4)
  )
  mu <- x@par$mu
  sigma <- x@par$sigma
  for (k in 1:K) {
    plot(mu[, k],
      type = "l", axes = F,
      col = "gray20", xlab = "", ylab = ""
    )
    axis(2, las = 2, cex.axis = .7)
    mtext(
      side = 2, las = 2, bquote(mu[k = .(k)]),
      cex = .6, line = 3
    )
  }
  for (k in 1:K) {
    plot(sigma[, k],
      type = "l", axes = F,
      col = "gray30", xlab = "", ylab = ""
    )
    axis(2, las = 2, cex.axis = .7)
    mtext(
      side = 2, las = 2, bquote(sigma[k = .(k)]),
      cex = .6, line = 3
    )
  }
  C <- x@hyper$C
  plot(c,
    type = "l", axes = F,
    col = "gray68", xlab = "", ylab = ""
  )
  axis(2, las = 2, cex.axis = 0.7)
  mtext(side = 2, las = 2, "C", cex = .6, line = 3)
  axis(1)
  mtext(side = 1, "Iterations", cex = .7, line = 3)
}

### --------------------------------------------------------------------
### .traces.Student.Hier
### @description    Plots traces for parameters of a univariate Student
###                 mixture.
### @par    x       an object of class mcmcoutputfix
###         dev     an object of class 'logical'
### @detail         Plots the traces for each component parameter of an
###                 Student mixture. If 'dev' is set to FALSE
###                 (TRUE is default) no device is created, instead
###                 the graphic can be stored to a file.
### @see            ?mcmcoutput, ?plotTraces
### @author         Lars Simon Zehnder
### --------------------------------------------------------------------
".traces.Student.Hier" <- function(x, dev) {
  K <- x@model@K
  trace.n <- 3 * K + 1
  if (.check.grDevice() && dev) {
    dev.new(title = "Traceplots")
  }
  par(
    mfrow = c(trace.n, 1), mar = c(1, 0, 0, 0),
    oma = c(4, 5, 4, 4)
  )
  mu <- x@par$mu
  sigma <- x@par$sigma
  df <- x@par$df
  for (k in 1:K) {
    plot(mu[, k],
      type = "l", axes = F,
      col = "gray20", xlab = "", ylab = ""
    )
    axis(2, las = 2, cex.axis = .7)
    mtext(
      side = 2, las = 2, bquote(mu[k = .(k)]),
      cex = .6, line = 3
    )
  }
  for (k in 1:K) {
    plot(sigma[, k],
      type = "l", axes = F,
      col = "gray30", xlab = "", ylab = ""
    )
    axis(2, las = 2, cex.axis = .7)
    mtext(
      side = 2, las = 2, bquote(sigma[k = .(k)]),
      cex = .6, line = 3
    )
  }
  for (k in 1:K) {
    plot(df[, k],
      type = "l", axes = F,
      col = "gray40", xlab = "", ylab = ""
    )
    axis(2, las = 2, cex.axis = .7)
    mtext(
      side = 2, las = 2, bquote(nu[k = .(k)]),
      cex = .6, line = 3
    )
  }
  C <- x@hyper$C
  plot(C,
    type = "l", axes = F,
    col = "gray68", xlab = "", ylab = ""
  )
  axis(2, las = 2, cex.axis = .7)
  mtext(
    side = 2, las = 2, "C", cex = .6,
    line = 3
  )
  axis(1)
  mtext(side = 1, "Iterations", cex = .7, line = 3)
}

"traces.Normult.Hier" <- function(x, dev, col) {
  .traces.Normult(x, dev, col)
  r <- x@model@r
  K <- x@model@K
  C <- x@hyper$C
  C.trace <- sapply(
    seq(1, x@M),
    function(i) sum(diag(qinmatr(C[i, ])))
  )
  C.logdet <- sapply(
    seq(1, x@M),
    function(i) log(det(qinmatr(C[i, ])))
  )
  # C traces
  mmax <- max(C.trace)
  mmin <- min(C.trace)
  if (.check.grDevice() && dev) {
    dev.new(title = "Traceplots Hyperparameters")
  }
  par(
    mfrow = c(2, 1), mar = c(1, 2, 0, 0),
    oma = c(4, 5, 4, 4)
  )
  if (col) {
    cscale <- rainbow(K, start = 0.5, end = 0)
  } else {
    cscale <- gray.colors(K, start = 0.5, end = 0.15)
  }
  plot(C.trace,
    type = "l", axes = F,
    col = cscale[K], xlab = "", ylab = ""
  )
  axis(2, las = 2, cex.axis = .7)
  mtext(
    side = 2, las = 2, bquote(tr(C)),
    cex = .6, line = 3
  )
  plot(C.logdet,
    type = "l", axes = F,
    col = cscale[K], xlab = "", ylab = ""
  )
  axis(2, las = 2, cex.axis = .7)
  name <- vector("character", K)
  mtext(
    side = 2, las = 2, bquote(log(det(C))),
    cex = .6, line = 3
  )
  axis(1)
  mtext(side = 1, "Iterations", cex = .7, line = 3)
}

".traces.Studmult.Hier" <- function(x, dev, col) {
  .traces.Studmult(x, dev, col)
  r <- x@model@r
  K <- x@model@K
  C <- x@hyper$C
  C.trace <- sapply(
    seq(1, x@M),
    function(i) sum(diag(qinmatr(C[i, ])))
  )
  C.logdet <- sapply(
    seq(1, x@M),
    function(i) log(det(qinmatr(C[i, ])))
  )

  # C traces
  mmax <- max(C.trace)
  mmin <- min(C.trace)
  if (.check.grDevice() && dev) {
    dev.new(title = "Traceplots Hyperparameters")
  }
  par(
    mfrow = c(2, 1), mar = c(1, 2, 0, 0),
    oma = c(4, 5, 4, 4)
  )
  if (col) {
    cscale <- rainbow(K, start = 0.5, end = 0)
  } else {
    cscale <- gray.colors(K, start = 0.5, end = 0.15)
  }
  plot(C.trace,
    type = "l", axes = F,
    col = cscale[K], xlab = "", ylab = ""
  )
  axis(2, las = 2, cex.axis = .7)
  mtext(
    side = 2, las = 2, bquote(tr(C)),
    cex = .6, line = 3
  )
  plot(C.logdet,
    type = "l", axes = F,
    col = cscale[K], xlab = "", ylab = ""
  )
  axis(2, las = 2, cex.axis = .7)
  mtext(
    side = 2, las = 2, bquote(log(det(C))),
    cex = .6, line = 3
  )
  axis(1)
  mtext(side = 1, "Iterations", cex = .7, line = 3)
}

### Plot Histograms
### Plot hist Poisson: Plots histograms for each component
### parameter and the hyper parameter 'b'.
".hist.Poisson.Hier" <- function(x, dev) {
  K <- x@model@K
  if (.check.grDevice() && dev) {
    dev.new(title = "Histograms")
  }
  lambda <- x@par$lambda
  b <- x@hyper$b
  vars <- cbind(lambda, b)
  if (K == 1) {
    lab.names <- list(bquote(lambda), "b")
    .symmetric.Hist(vars, lab.names)
  } else {
    lab.names <- vector("list", K + 1)
    for (k in 1:K) {
      lab.names[[k]] <- bquote(lambda[.(k)])
    }
    lab.names[[K + 1]] <- "b"
    .symmetric.Hist(vars, lab.names)
  }
}

".hist.Normal.Hier" <- function(x, dev) {
  K <- x@model@K
  mu <- x@par$mu
  sigma <- x@par$sigma
  C <- x@hyper$C
  if (K == 1) {
    if (.check.grDevice() && dev) {
      dev.new(title = "Histogram Mu")
    }
    .symmetric.Hist(mu, list(bquote(mu)))
    if (.check.grDevice() && dev) {
      dev.new(title = "Histogram Sigma")
    }
    .symmetric.Hist(sigma, list(bquote(sigma)))
  } else {
    mu.lab.names <- vector("list", K)
    sigma.lab.names <- vector("list", K)
    for (k in 1:K) {
      mu.lab.names[[k]] <- bquote(mu[.(k)])
      sigma.lab.names[[k]] <- bquote(sigma[.(k)])
    }
    if (.check.grDevice() && dev) {
      dev.new(title = "Histograms Mu")
    }
    .symmetric.Hist(mu, mu.lab.names)
    if (.check.grDevice() && dev) {
      dev.new(title = "Histograms Sigma")
    }
    .symmetric.Hist(sigma, sigma.lab.names)
  }
  if (.check.grDevice() && dev) {
    dev.new(title = "Histogram Hyperparameter C")
  }
  .symmetric.Hist(C, "C")
}

".hist.Student.Hier" <- function(x, dev) {
  K <- x@model@K
  mu <- x@par$mu
  sigma <- x@par$sigma
  degf <- x@par$df
  C <- x@hyper$C
  if (K == 1) {
    if (.check.grDevice() && dev) {
      dev.new(title = "Histogram Mu")
    }
    .symmetric.Hist(mu, list(bquote(mu)))
    if (.check.grDevice() && dev) {
      dev.new(title = "Histogram Sigma")
    }
    .symmetric.Hist(sigma, list(bquote(sigma)))
    if (.check.grDevice() && dev) {
      dev.new(title = "Histogram Degrees of Freedom")
    }
    .symmetric.Hist(degf, list(bquote(nu)))
  } else {
    mu.lab.names <- vector("list", K)
    sigma.lab.names <- vector("list", K)
    degf.lab.names <- vector("list", K)
    for (k in 1:K) {
      mu.lab.names[[k]] <- bquote(mu[.(k)])
      sigma.lab.names[[k]] <- bquote(sigma[.(k)])
      degf.lab.names[[k]] <- bquote(nu[.(k)])
    }
    if (.check.grDevice() && dev) {
      dev.new(title = "Histograms Mu")
    }
    .symmetric.Hist(mu, mu.lab.names)
    if (.check.grDevice() && dev) {
      dev.new(title = "Histograms Sigma")
    }
    .symmetric.Hist(sigma, sigma.lab.names)
    if (.check.grDevice() && dev) {
      dev.new(title = "Histograms Degrees of Freedom")
    }
    .symmetric.Hist(degf, degf.lab.names)
  }
  if (.check.grDevice() && dev) {
    dev.new(title = "Histogram Hyperparameter C")
  }
  .symmetric.Hist(C, "C")
}

".hist.Normult.Hier" <- function(x, dev) {
  K <- x@model@K
  r <- x@model@r
  mu <- x@par$mu
  sigma <- x@par$sigma
  logdetC <- sapply(seq(1, x@M), function(i) log(det(qinmatr(x@hyper$C[i, ]))))
  trC <- sapply(seq(1, x@M), function(i) sum(diag(qinmatr(x@hyper$C[i, ]))))
  for (rr in 1:r) {
    if (K == 1) {
      if (.check.grDevice() & dev) {
        dev.new(title = paste("Histograms Feature ", rr,
          " Mu",
          sep = ""
        ))
      }
      .symmetric.Hist(mu[, rr, ], list(bquote(mu)))
      if (.check.grDevice() & dev) {
        dev.new(title = paste("Histograms Feature ", rr,
          " Sigma",
          sep = ""
        ))
      }
      .symmetric.Hist(sigma[, rr, ], list(bquote(sigma)))
    } else {
      mu.lab.names <- vector("list", K)
      sigma.lab.names <- vector("list", K)
      for (k in 1:K) {
        mu.lab.names[[k]] <- bquote(mu[.(k)])
        sigma.lab.names[[k]] <- bquote(sigma[.(k)])
      }
      if (.check.grDevice() & dev) {
        dev.new(title = paste("Histograms Feature ", rr,
          " Mu",
          sep = ""
        ))
      }
      .symmetric.Hist(mu[, rr, ], mu.lab.names)
      if (.check.grDevice() & dev) {
        dev.new(title = paste("Histograms Feature ", rr,
          " Sigma",
          sep = ""
        ))
      }
      .symmetric.Hist(sigma[, rr, ], sigma.lab.names)
    }
  }
  if (.check.grDevice() && dev) {
    dev.new(title = "Histograms Hyperparameter C")
  }
  C.lab.names <- vector("list", 2)
  C.lab.names[[1]] <- "log(det(C))"
  C.lab.names[[2]] <- "tr(C)"
  .symmetric.Hist(cbind(logdetC, trC), C.lab.names)
}

".hist.Studmult.Hier" <- function(x, dev) {
  K <- x@model@K
  r <- x@model@r
  mu <- x@par$mu
  sigma <- x@par$sigma
  degf <- x@par$df
  logdetC <- sapply(seq(1, x@M), function(i) log(det(qinmatr(x@hyper$C[i, ]))))
  trC <- sapply(seq(1, x@M), function(i) sum(diag(qinmatr(x@hyper$C[i, ]))))
  for (rr in 1:r) {
    if (K == 1) {
      if (.check.grDevice() & dev) {
        dev.new(title = paste("Histograms Feature ", rr,
          " Mu",
          sep = ""
        ))
      }
      .symmetric.Hist(mu[, rr, ], list(bquote(mu)))
      if (.check.grDevice() & dev) {
        dev.new(title = paste("Histograms Feature ", rr,
          " Sigma",
          sep = ""
        ))
      }
      .symmetric.Hist(sigma[, rr, ], list(bquote(sigma)))
    } else {
      mu.lab.names <- vector("list", K)
      sigma.lab.names <- vector("list", K)
      for (k in 1:K) {
        mu.lab.names[[k]] <- bquote(mu[.(k)])
        sigma.lab.names[[k]] <- bquote(sigma[.(k)])
      }
      if (.check.grDevice() & dev) {
        dev.new(title = paste("Histograms Feature ", rr,
          " Mu",
          sep = ""
        ))
      }
      .symmetric.Hist(mu[, rr, ], mu.lab.names)
      if (.check.grDevice() & dev) {
        dev.new(title = paste("Histograms Feature ", rr,
          " Sigma",
          sep = ""
        ))
      }
      .symmetric.Hist(sigma[, rr, ], sigma.lab.names)
    }
  }
  if (K == 1) {
    if (.check.grDevice() & dev) {
      dev.new(title = paste("Histograms Feature ", rr,
        " Mu",
        sep = ""
      ))
    }
    .symmetric.Hist(degf[, rr, ], list(bquote(nu)))
  } else {
    degf.lab.names <- vector("list", K)
    for (k in 1:K) {
      degf.lab.names[[k]] <- bquote(nu[.(k)])
    }
    if (.check.grDevice() & dev) {
      dev.new(title = paste("Histograms Feature ", rr,
        " Sigma",
        sep = ""
      ))
    }
    .symmetric.Hist(degf[, rr, ], degf.lab.names)
  }
  if (.check.grDevice() && dev) {
    dev.new(title = "Histograms Hyperparameter C")
  }
  C.lab.names <- vector("list", 2)
  C.lab.names[[1]] <- "log(det(C))"
  C.lab.names[[2]] <- "tr(C)"
  .symmetric.Hist(cbind(logdetC, trC), C.lab.names)
}

### Plot Densities
### Plot Dens Poisson Hier: Plots Kernel densities for each
### component parameter and the hyper parameter 'b'.
".dens.Poisson.Hier" <- function(x, dev) {
  K <- x@model@K
  if (.check.grDevice() && dev) {
    dev.new("Densities")
  }
  lambda <- x@par$lambda
  b <- x@hyper$b
  vars <- cbind(lambda, b)
  if (K == 1) {
    lab.names <- list(bquote(lambda), "b")
    .symmetric.Dens(vars, lab.names)
  } else {
    lab.names <- vector("list", K + 1)
    for (k in seq(1, K)) {
      lab.names[[k]] <- bquote(lambda[.(k)])
    }
    lab.names[[K + 1]] <- "b"
    .symmetric.Dens(vars, lab.names)
  }
}

".dens.Normal.Hier" <- function(x, dev) {
  K <- x@model@K
  mu <- x@par$mu
  sigma <- x@par$sigma
  C <- x@hyper$C
  if (K == 1) {
    if (.check.grDevice() && dev) {
      dev.new(title = "Density Mu")
    }
    .symmetric.Dens(mu, list(bquote(mu)))
    if (.check.grDevice() && dev) {
      dev.new(title = "Density Sigma")
    }
    .symmetric.Dens(sigma, list(bquote(sigma)))
  } else {
    mu.lab.names <- vector("list", K)
    sigma.lab.names <- vector("list", K)
    for (k in 1:K) {
      mu.lab.names[[k]] <- bquote(mu[.(k)])
      sigma.lab.names[[k]] <- bquote(sigma[.(k)])
    }
    if (.check.grDevice() && dev) {
      dev.new(title = "Densities Mu")
    }
    .symmetric.Dens(mu, mu.lab.names)
    if (.check.grDevice() && dev) {
      dev.new(title = "Densities Sigma")
    }
    .symmetric.Dens(sigma, sigma.lab.names)
  }
  if (.check.grDevice() && dev) {
    dev.new(title = "Histogram Hyperparameter C")
  }
  .symmetric.Dens(C, "C")
}

".dens.Student.Hier" <- function(x, dev) {
  K <- x@model@K
  mu <- x@par$mu
  sigma <- x@par$sigma
  degf <- x@par$df
  C <- x@hyper$C
  if (K == 1) {
    if (.check.grDevice() && dev) {
      dev.new(title = "Density Mu")
    }
    .symmetric.Dens(mu, list(bquote(mu)))
    if (.check.grDevice() && dev) {
      dev.new(title = "Density Sigma")
    }
    .symmetric.Dens(sigma, list(bquote(sigma)))
    if (.check.grDevice() && dev) {
      dev.new(title = "Density Degrees of Freedom")
    }
    .symmetric.Dens(degf, list(bquote(nu)))
  } else {
    mu.lab.names <- vector("list", K)
    sigma.lab.names <- vector("list", K)
    degf.lab.names <- vector("list", K)
    for (k in 1:K) {
      mu.lab.names[[k]] <- bquote(mu[.(k)])
      sigma.lab.names[[k]] <- bquote(sigma[.(k)])
      degf.lab.names[[k]] <- bquote(nu[.(k)])
    }
    if (.check.grDevice() && dev) {
      dev.new(title = "Densities Mu")
    }
    .symmetric.Dens(mu, mu.lab.names)
    if (.check.grDevice() && dev) {
      dev.new(title = "Densities Sigma")
    }
    .symmetric.Dens(sigma, sigma.lab.names)
    if (.check.grDevice() && dev) {
      dev.new(title = "Densities Degrees of Freedom")
    }
    .symmetric.Dens(degf, degf.lab.names)
  }
  if (.check.grDevice() && dev) {
    dev.new(title = "Density Hyperparameter C")
  }
  .symmetric.Dens(C, "C")
}

"dens.Normult.Hier" <- function(x, dev) {
  K <- x@model@K
  r <- x@model@r
  mu <- x@par$mu
  sigma <- x@par$sigma
  logdetC <- sapply(seq(1, x@M), function(i) log(det(qinmatr(x@hyper$C[i, ]))))
  trC <- sapply(seq(1, x@M), function(i) sum(diag(qinmatr(x@hyper$C[i, ]))))
  for (rr in 1:r) {
    if (K == 1) {
      if (.check.grDevice() & dev) {
        dev.new(title = paste("Densities Feature ", rr,
          " Mu",
          sep = ""
        ))
      }
      .symmetric.Dens(mu[, rr, ], list(bquote(mu)))
      if (.check.grDevice() & dev) {
        dev.new(title = paste("Densities Feature ", rr,
          " Sigma",
          sep = ""
        ))
      }
      .symmetric.Dens(sigma[, rr, ], list(bquote(sigma)))
    } else {
      mu.lab.names <- vector("list", K)
      sigma.lab.names <- vector("list", K)
      for (k in 1:K) {
        mu.lab.names[[k]] <- bquote(mu[.(k)])
        sigma.lab.names[[k]] <- bquote(sigma[.(k)])
      }
      if (.check.grDevice() & dev) {
        dev.new(title = paste("Densities Feature ", rr,
          " Mu",
          sep = ""
        ))
      }
      .symmetric.Dens(mu[, rr, ], mu.lab.names)
      if (.check.grDevice() & dev) {
        dev.new(title = paste("Densities Feature ", rr,
          " Sigma",
          sep = ""
        ))
      }
      .symmetric.Dens(sigma[, rr, ], sigma.lab.names)
    }
  }
  if (.check.grDevice() && dev) {
    dev.new(title = "Densities Hyperparameter C")
  }
  C.lab.names <- vector("list", 2)
  C.lab.names[[1]] <- "log(det(C))"
  C.lab.names[[2]] <- "tr(C)"
  .symmetric.Dens(cbind(logdetC, trC), C.lab.names)
}

"dens.Studmult.Hier" <- function(x, dev) {
  K <- x@model@K
  r <- x@model@r
  mu <- x@par$mu
  sigma <- x@par$sigma
  degf <- x@par$df
  logdetC <- sapply(seq(1, x@M), function(i) log(det(qinmatr(x@hyper$C[i, ]))))
  trC <- sapply(seq(1, x@M), function(i) sum(diag(qinmatr(x@hyper$C[i, ]))))
  for (rr in 1:r) {
    if (K == 1) {
      if (.check.grDevice() & dev) {
        dev.new(title = paste("Densities Feature ", rr,
          " Mu",
          sep = ""
        ))
      }
      .symmetric.Dens(mu[, rr, ], list(bquote(mu)))
      if (.check.grDevice() & dev) {
        dev.new(title = paste("Densities Feature ", rr,
          " Sigma",
          sep = ""
        ))
      }
      .symmetric.Dens(sigma[, rr, ], list(bquote(sigma)))
    } else {
      mu.lab.names <- vector("list", K)
      sigma.lab.names <- vector("list", K)
      for (k in 1:K) {
        mu.lab.names[[k]] <- bquote(mu[.(k)])
        sigma.lab.names[[k]] <- bquote(sigma[.(k)])
      }
      if (.check.grDevice() & dev) {
        dev.new(title = paste("Densities Feature ", rr,
          " Mu",
          sep = ""
        ))
      }
      .symmetric.Dens(mu[, rr, ], mu.lab.names)
      if (.check.grDevice() & dev) {
        dev.new(title = paste("Densities Feature ", rr,
          " Sigma",
          sep = ""
        ))
      }
      .symmetric.Dens(sigma[, rr, ], sigma.lab.names)
    }
  }
  if (K == 1) {
    if (.check.grDevice() & dev) {
      dev.new(title = "Density Degrees of Freedom")
    }
    .symmetric.Dens(degf[, rr, ], list(bquote(nu)))
  } else {
    degf.lab.names <- vector("list", K)
    for (k in 1:K) {
      degf.lab.names[[k]] <- bquote(nu[.(k)])
    }
    if (.check.grDevice() & dev) {
      dev.new(title = "Densities Degrees of Freedom")
    }
    .symmetric.Dens(degf[, rr, ], degf.lab.names)
  }
  if (.check.grDevice() && dev) {
    dev.new(title = "Densities Hyperparameter C")
  }
  C.lab.names <- vector("list", 2)
  C.lab.names[[1]] <- "log(det(C))"
  C.lab.names[[2]] <- "tr(C)"
  .symmetric.Dens(cbind(logdetC, trC), C.lab.names)
}

### Logic
### Logic subseq Hier: Creates a subsequence for the sample
### of the Poisson hyper parameter 'b'.
".subseq.Poisson.Hier" <- function(obj, index) {
  obj@hyper$b <- array(obj@hyper$b[index],
    dim = c(obj@M, 1)
  )
  return(obj)
}

".subseq.Norstud.Hier" <- function(obj, index) {
  obj@hyper$C <- array(obj@hyper$C[index],
    dim = c(obj@M, 1)
  )
  return(obj)
}

".subseq.Normultstud.Hier" <- function(obj, index) {
  obj@hyper$C <- array(obj@hyper$C[index, ],
    dim = c(obj@M, obj@model@K)
  )
  return(obj)
}
