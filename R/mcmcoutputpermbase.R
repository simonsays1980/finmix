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

.mcmcoutputpermbase <- setClass("mcmcoutputpermbase",
  contains = c(
    "mcmcpermind",
    "mcmcoutputbase"
  ),
  validity = function(object) {
    ## else: OK
    TRUE
  }
)

setMethod(
  "initialize", "mcmcoutputpermbase",
  function(.Object, mcmcoutput, Mperm = integer(),
           parperm = list(), relabel = character(),
           weightperm = array(), logperm = list(),
           entropyperm = array(), STperm = array(),
           Sperm = array(), NKperm = array()) {
    .Object@M <- mcmcoutput@M
    .Object@burnin <- mcmcout@burnin
    .Object@ranperm <- mcmcoutput@ranperm
    .Object@par <- mcmcoutput@par
    .Object@weight <- mcmcoutput@weight
    .Object@log <- mcmcoutput@log
    .Object@ST <- mcmcoutput@ST
    .Object@S <- mcmcoutput@S
    .Object@NK <- mcmcoutput@NK
    .Object@clust <- mcmcoutput@clust
    .Object@model <- mcmcoutput@model
    .Object@prior <- mcmcoutput@prior
    .Object@Mperm <- Mperm
    .Object@parperm <- parperm
    .Object@relabel <- relabel
    .Object@weightperm <- weightperm
    .Object@logperm <- logperm
    .Object@entropyperm <- entropyperm
    .Object@STperm <- STperm
    .Object@Sperm <- Sperm
    .Object@NKperm <- NKperm
    .Object
  }
)

setMethod(
  "show", "mcmcoutputpermbase",
  function(object) {
    cat("Object 'mcmcoutputperm'\n")
    cat("     class       :", class(object), "\n")
    cat("     M           :", object@M, "\n")
    cat("     burnin      :", object@burnin, "\n")
    cat("     ranperm     :", object@ranperm, "\n")
    cat("     relabel     :", object@relabel, "\n")
    cat(
      "     par         : List of",
      length(object@par), "\n"
    )
    cat(
      "     log         : List of",
      length(object@log), "\n"
    )
    cat(
      "     ST          :",
      paste(dim(object@ST), collapse = "x"), "\n"
    )
    if (!all(is.na(object@S))) {
      cat(
        "     S           :",
        paste(dim(object@S), collapse = "x"), "\n"
      )
    }
    cat(
      "     NK          :",
      paste(dim(object@NK), collapse = "x"), "\n"
    )
    cat(
      "     clust       :",
      paste(dim(object@clust), collapse = "x"), "\n"
    )
    cat("     Mperm       :", object@Mperm, "\n")
    cat(
      "     parperm     : List of",
      length(object@parperm), "\n"
    )
    cat(
      "     weightperm  :",
      paste(dim(object@weightperm), collapse = "x"), "\n"
    )
    cat(
      "     logperm     : List of",
      length(object@logperm), "\n"
    )
    cat(
      "     entropyperm :",
      paste(dim(object@entropyperm), collapse = "x"), "\n"
    )
    cat(
      "     STperm      :",
      paste(dim(object@STperm), collapse = "x"), "\n"
    )
    if (!all(is.na(object@Sperm))) {
      cat(
        "     Sperm       :",
        paste(dim(object@Sperm), collapse = "x"), "\n"
      )
    }
    cat(
      "     NKperm      :",
      paste(dim(object@NKperm), collapse = "x"), "\n"
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
    x = "mcmcoutputpermbase",
    dev = "ANY",
    lik = "ANY",
    col = "ANY"
  ),
  function(x, dev = TRUE, lik = 1, col = FALSE, ...) {
    dist <- x@model@dist
    if (lik %in% c(0, 1)) {
      if (dist == "poisson") {
        .permtraces.Poisson.Base(x, dev)
      } else if (dist == "binomial") {
        .permtraces.Binomial.Base(x, dev)
      } else if (dist == "exponential") {
        .permtraces.Exponential.Base(x, dev)
      } else if (dist == "normal") {
        .permtraces.Normal(x, dev)
        .permtraces.Weights.Base(x, dev, col)
      } else if (dist == "student") {
        .permtraces.Student(x, dev)
        .permtraces.Weights.Base(x, dev, col)
      } else if (dist == "normult") {
        .permtraces.Normult(x, dev, col)
        .permtraces.Weights.Base(x, dev, col)
      } else if (dist == "studmult") {
        .permtraces.Studmult(x, dev, col)
        .permtraces.Weights.Base(x, dev, col)
      }
    }
    if (lik %in% c(1, 2)) {
      ## log ##
      .permtraces.Log.Base(x, dev)
    }
  }
)

setMethod(
  "plotHist", signature(
    x = "mcmcoutputpermbase",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    dist <- x@model@dist
    if (dist == "poisson") {
      .permhist.Poisson.Base(x, dev)
    } else if (dist == "binomial") {
      .permhist.Binomial.Base(x, dev)
    }
  }
)

setMethod(
  "plotDens", signature(
    x = "mcmcoutputpermbase",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    dist <- x@model@dist
    if (dist == "poisson") {
      .permdens.Poisson.Base(x, dev)
    } else if (dist == "binomial") {
      .permdens.Binomial.Base(x, dev)
    }
  }
)

setMethod(
  "plotPointProc", signature(
    x = "mcmcoutputpermbase",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    dist <- x@model@dist
    if (dist == "poisson") {
      .permpointproc.Poisson(x, dev)
    } else if (dist == "binomial") {
      .permpointproc.Binomial(x, dev)
    }
  }
)

setMethod(
  "plotSampRep", signature(
    x = "mcmcoutputpermbase",
    dev = "ANY"
  ),
  function(x, dev, ...) {
    dist <- x@model@dist
    if (dist == "poisson") {
      .permsamprep.Poisson(x, dev)
    } else if (dist == "binomial") {
      .permsamprep.Binomial(x, dev)
    }
  }
)

setMethod(
  "plotPostDens", signature(
    x = "mcmcoutputpermbase",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    dist <- x@model@dist
    if (dist == "poisson") {
      .permpostdens.Poisson(x, dev)
    } else if (dist == "binomial") {
      .permpostdens.Binomial(x, dev)
    }
  }
)

### Private functions.
### These functions are not exported.

### Plot
### Traces
### Traces Poisson: Plots traces for all Poisson parameters
### and the weights.
".permtraces.Poisson.Base" <- function(x, dev) {
  K <- x@model@K
  trace.n <- K * 2 - 1
  if (.check.grDevice() && dev) {
    dev.new(title = "Traceplots (permuted)")
  }
  par(
    mfrow = c(trace.n, 1), mar = c(1, 0, 0, 0),
    oma = c(4, 5, 4, 4)
  )
  lambda <- x@parperm$lambda
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
  weight <- x@weightperm
  for (k in 1:(K - 1)) {
    plot(weight[, k],
      type = "l", axes = F,
      col = "gray47", xlab = "", ylab = ""
    )
    axis(2, las = 2, cex.axis = 0.7)
    mtext(
      side = 2, las = 2, bquote(eta[k = .(k)]),
      cex = 0.6, line = 3
    )
  }
  axis(1)
  mtext(side = 1, "Iterations", cex = 0.7, line = 3)
}

".permtraces.Binomial.Base" <- function(x, dev) {
  K <- x@model@K
  trace.n <- K * 2 - 1
  if (.check.grDevice() && dev) {
    dev.new(title = "Traceplots (permuted)")
  }
  par(
    mfrow = c(trace.n, 1), mar = c(1, 0, 0, 0),
    oma = c(4, 5, 4, 4)
  )
  p <- x@parperm$p
  for (k in 1:K) {
    plot(p[, k],
      type = "l", axes = F,
      col = "gray20", xlab = "", ylab = ""
    )
    axis(2, las = 2, cex.axis = 0.7)
    mtext(
      side = 2, las = 2, bquote(p[k = .(k)]),
      cex = 0.6, line = 3
    )
  }
  weight <- x@weightperm
  for (k in 1:(K - 1)) {
    plot(weight[, k],
      type = "l", axes = F,
      col = "gray47", xlab = "", ylab = ""
    )
    axis(2, las = 2, cex.axis = 0.7)
    mtext(
      side = 2, las = 2, bquote(eta[k = .(k)]),
      cex = 0.6, line = 3
    )
  }
  axis(1)
  mtext(side = 1, "Iterations", cex = 0.7, line = 3)
}

".permtraces.Exponential.Base" <- function(x, dev) {
  K <- x@model@K
  trace.n <- K * 2 - 1
  if (.check.grDevice() && dev) {
    dev.new(title = "Traceplots")
  }
  par(
    mfrow = c(trace.n, 1), mar = c(1, 0, 0, 0),
    oma = c(4, 5, 4, 4)
  )
  lambda <- x@parperm$lambda
  for (k in 1:K) {
    plot(lambda[, k],
      type = "l", axes = F,
      col = "gray20", xlab = "", ylab = ""
    )
    axis(2, las = 2, cex.axis = .7)
    mtext(
      side = 2, las = 2, bquote(lambda[k = .(k)]),
      cex = .6, line = 3
    )
  }
  weight <- x@weight
  for (k in 1:(K - 1)) {
    plot(weight[, k],
      type = "l", axes = F,
      col = "gray47", xlab = "", ylab = ""
    )
    axis(2, las = 2, cex.axis = .7)
    mtext(
      side = 2, las = 2, bquote(eta[k = .(k)]),
      cex = .6, line = 3
    )
  }
  axis(1)
  mtext(side = 1, "Iterations", cex = .7, line = 3)
}

".permtraces.Weights.Base" <- function(x, dev, col) {
  weight <- x@weightperm
  K <- x@model@K
  if (.check.grDevice() && dev) {
    dev.new(title = "Traceplots Weights")
  }
  if (col) {
    cscale <- rainbow(K, start = 0.5, end = 0)
  } else {
    cscale <- gray.colors(K, start = 0.5, end = 0.15)
  }

  plot(weight[, 1],
    type = "l", axes = F,
    col = cscale[1], xlab = "", ylab = "", ylim = c(0, 1.2)
  )
  for (k in 2:K) {
    lines(weight[, k], col = cscale[k])
  }
  axis(2, las = 2, cex.axis = .7)
  mtext(
    side = 2, las = 2, bquote(eta),
    cex = .6, line = 3
  )
  name <- vector("character", K)
  for (k in 1:K) {
    name[k] <- paste("k = ", k, sep = "")
  }
  legend("top",
    legend = name, col = cscale, lty = 1,
    horiz = TRUE, cex = .7
  )
  axis(1)
  mtext(side = 1, "Iterations", cex = .7, line = 3)
}

### Traces log-likelihoods: Plots traces for the log-likelihoods.
".permtraces.Log.Base" <- function(x, dev) {
  if (.check.grDevice() && dev) {
    dev.new(title = "Log Likelihood Traceplots (permuted)")
  }
  if (col) {
    cscale <- rainbow(3, start = 0, end = .5)
  } else {
    cscale <- gray.colors(3, start = 0, end = .15)
  }
  par(
    mfrow = c(3, 1), mar = c(1, 0, 0, 0),
    oma = c(4, 5, 4, 4)
  )
  mixlik <- x@logperm$mixlik
  plot(mixlik,
    type = "l", axes = F,
    col = cscale[3], xlab = "", ylab = ""
  )
  axis(2, las = 2, cex.axis = 0.7)
  mtext(
    side = 2, las = 3, "mixlik", cex = 0.6,
    line = 3
  )
  mixprior <- x@logperm$mixprior
  plot(mixprior,
    type = "l", axes = F,
    col = cscale[2], xlab = "", ylab = ""
  )
  axis(2, las = 2, cex.axis = 0.7)
  mtext(
    side = 2, las = 3, "mixprior", cex = 0.6,
    line = 3
  )
  cdpost <- x@logperm$cdpost
  plot(mixprior,
    type = "l", axes = F,
    col = cscale[3], xlab = "", ylab = ""
  )
  axis(2, las = 2, cex.axis = 0.7)
  mtext(
    side = 2, las = 3, "cdpost", cex = 0.6,
    line = 3
  )
  axis(1)
  mtext(side = 1, "Iterations", cex = 0.7, line = 3)
}

### Histograms
### Histograms Poisson: Plots histograms for all Poisson
### parameters and the weights.
".permhist.Poisson.Base" <- function(x, dev) {
  K <- x@model@K
  if (.check.grDevice() && dev) {
    dev.new(title = "Histograms (permuted)")
  }
  lambda <- x@parperm$lambda
  weight <- x@weightperm
  vars <- cbind(lambda, weight[, seq(1:(K - 1))])
  lab.names <- vector("list", 2 * K - 1)
  for (k in 1:K) {
    lab.names[[k]] <- bquote(lambda[.(k)])
  }
  for (k in (K + 1):(2 * K - 1)) {
    lab.names[[k]] <- bquote(eta[.(k - K)])
  }
  .symmetric.Hist(vars, lab.names)
}

".permhist.Binomial.Base" <- function(x, dev) {
  K <- x@model@K
  if (.check.grDevice() && dev) {
    dev.new(title = "Histograms (permuted)")
  }
  p <- x@parperm$p
  weight <- x@weightperm
  vars <- cbind(p, weight[, seq(1:(K - 1))])
  lab.names <- vector("list", 2 * K - 1)
  for (k in 1:K) {
    lab.names[[k]] <- bquote(p[.(k)])
  }
  for (k in (K + 1):(2 * K - 1)) {
    lab.names[[k]] <- bquote(eta[.(k - K)])
  }
  .symmetric.Hist(vars, lab.names)
}

### Densities
### Densities Poisson: Plots Kernel densities for all Poisson
### parameters and weights.
".permdens.Poisson.Base" <- function(x, dev) {
  K <- x@model@K
  if (.check.grDevice() && dev) {
    dev.new(title = "Histograms (permuted)")
  }
  lambda <- x@parperm$lambda
  weight <- x@weightperm
  vars <- cbind(lambda, weight[, seq(1:(K - 1))])
  lab.names <- vector("list", 2 * K - 1)
  for (k in 1:K) {
    lab.names[[k]] <- bquote(lambda[.(k)])
  }
  for (k in (K + 1):(2 * K - 1)) {
    lab.names[[k]] <- bquote(eta[.(k - K)])
  }
  .symmetric.Dens(vars, lab.names)
}

".permdens.Binomial.Base" <- function(x, dev) {
  K <- x@model@K
  if (.check.grDevice() && dev) {
    dev.new(title = "Histograms (permuted)")
  }
  p <- x@parperm$p
  weight <- x@weightperm
  vars <- cbind(p, weight[, seq(1:(K - 1))])
  lab.names <- vector("list", 2 * K - 1)
  for (k in 1:K) {
    lab.names[[k]] <- bquote(p[.(k)])
  }
  for (k in (K + 1):(2 * K - 1)) {
    lab.names[[k]] <- bquote(eta[.(k - K)])
  }
  .symmetric.Dens(vars, lab.names)
}
