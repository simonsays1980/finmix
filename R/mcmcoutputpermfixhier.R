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

#' Finmix `mcmcoutputpermfixhier` class 
#' 
#' @description 
#' This class defines objects to store the outputs from permuting the MCMC 
#' samples. Due to label switching the sampled component parameters are usually 
#' not assigned to the same component in each iteration. To overcome this issue 
#' the samples are permuted by using a relabeling algorithm (usually K-means) 
#' to reassign parameters. Note that due to assignment of parameters from the 
#' same iteration to the same component, the sample size could shrink. 
#' 
#' This class stores the permuted parameters together with the new MCMC sample 
#' size and the mixture log-likelihood, the prior log-likelihood, and the 
#' complete data posterior log-likelihood.
#' 
#' Note this class inherits all slots from its parent classes.  
#' 
#' @exportClass mcmcoutputpermfixhier
#' @rdname mcmcoutputpermfixhier-class
#' @seealso 
#' * [mcmcoutputpermfix-class] for the parent class
#' * [mcmcpermfix-class] for the parent class
#' * [mcmcpermute()] for performing permutation of MCMC samples
.mcmcoutputpermfixhier <- setClass("mcmcoutputpermfixhier",
  contains = c("mcmcpermfixhier", "mcmcoutputfixhier"),
  validity = function(object) {
    ## else: OK
    TRUE
  }
)

#' Initializer of the `mcmcoutputpermfixhier` class
#' 
#' @description
#' Only used implicitly. The initializer stores the data into the slots of the 
#' passed-in object.
#' 
#' @param .Object An object: see the "initialize Methods" section in 
#'   [initialize].
#' @param mcmcoutput An `mcmcoutput` class containing the results from MCMC 
#'   sampling.
#' @param Mperm An integer defining the number of permuted MCMC samples.
#' @param parperm A named list containing the permuted component parameter 
#'   samples from MCMC sampling
#' @param logperm A named list containing the mixture log-likelihood, the 
#'   prior log-likelihood, and the complete data posterior log-likelihood 
#'   for the permuted MCMC samples.
#' @param hyperperm A named list containing the permuted parameters of the 
#'   hierarchical prior.
#' 
#' @keywords internal
#' 
#' @seealso 
#' * [Classes_Details] for details of class definitions, and 
#' * [setOldClass] for the relation to S3 classes
setMethod(
  "initialize", "mcmcoutputpermfixhier",
  function(.Object, mcmcoutput, Mperm = integer(),
           parperm = list(), logperm = list(), hyperperm = list()) {
    .Object@M <- mcmcoutput@M
    .Object@burnin <- mcmcoutput@burnin
    .Object@ranperm <- mcmcoutput@ranperm
    .Object@par <- mcmcoutput@par
    .Object@log <- mcmcoutput@log
    .Object@hyper <- mcmcoutput@hyper
    .Object@model <- mcmcoutput@model
    .Object@prior <- mcmcoutput@prior
    .Object@Mperm <- Mperm
    .Object@parperm <- parperm
    .Object@logperm <- logperm
    .Object@hyperperm <- hyperperm
    .Object
  }
)

#' Shows a summary of an `mcmcoutputpermfixhier` object.
#' 
#' @description
#' Calling [show()] on an `mcmcoutputpermfixhier` object gives an overview 
#' of the `mcmcoutputpermfixhier` object.
#' 
#' @param object An `mcmcoutputpermfixhier` object.
#' @returns A console output listing the slots and summary information about
#'   each of them. 
#' @exportMethod show
#' @noRd
setMethod(
  "show", "mcmcoutputpermfixhier",
  function(object) {
    cat("Object 'mcmcoutputperm'\n")
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
      "     hyper       : List of",
      length(object@hyper), "\n"
    )
    cat("     Mperm       :", object@Mperm, "\n")
    cat(
      "     parperm     : List of",
      length(object@parperm), "\n"
    )
    cat(
      "     logperm     : List of",
      length(object@logperm), "\n"
    )
    cat(
      "     hyperperm   : List of",
      length(object@hyperperm), "\n"
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

#' Plot traces of MCMC sampling
#' 
#' @description 
#' Calling [plotTraces()] plots the MCMC traces of the mixture log-likelihood 
#' , the mixture log-likelihood of the prior distribution, the log-likelihood 
#' of the complete data posterior, or the weights and parameters if `lik` is 
#' set to `1`.s 
#' 
#' If `lik` is set to `0` the parameters of the components and the posterior 
#' parameters are plotted together with `K-1` weights.
#' 
#' @param x An `mcmcoutputpermfixhier` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param lik An integer indicating, if the log-likelihood traces should be 
#'   plotted (default). If set to `0` the traces for the parameters 
#'   and weights are plotted instead. 
#' @param col A logical indicating, if the plot should be colored.
#' @param ... Further arguments to be passed to the plotting function.
#' @return A plot of the traces of the MCMC samples.
#' @exportMethod plotTraces
#' 
#' @examples 
#' \dontrun{
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2, 
#'                  indicfix = TRUE)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc(storepost = FALSE)
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' f_outputperm <- mcmcpermute(f_output)
#' plotTraces(f_outputperm, lik = 0)
#' }
#' 
#' @seealso 
#' * [mixturemcmc()] for performing MCMC sampling
#' * [mcmcpermute()] for permuting MCMC samples
#' * [plotHist()] for plotting histograms of sampled values
#' * [plotDens()] for plotting densities of sampled values
#' * [plotSampRep()] for plotting sampling representations of sampled values
#' * [plotPointProc()] for plotting point processes for sampled values
#' * [plotPostDens()] for plotting the posterior density of component parameters
setMethod(
  "plotTraces", signature(
    x = "mcmcoutputpermfixhier",
    dev = "ANY",
    lik = "ANY",
    col = "ANY"
  ),
  function(x, dev = TRUE, lik = 1, col = FALSE, ...) {
    dist <- x@model@dist
    if (lik %in% c(0, 1)) {
      if (dist == "poisson") {
        .permtraces.Poisson.Hier(x, dev)
      } else if (dist == "binomial") {
        .permtraces.Binomial(x, dev)
      } else if (dist == "exponential") {
        callNextMethod(x, dev, lik, col, ...)
      } else if (dist == "normal") {
        .permtraces.Normal.Hier(x, dev)
      } else if (dist == "student") {
        .permtraces.Student.Hier(x, dev)
      } else if (dist == "normult") {
        .permtraces.Normult.Hier(x, dev, col)
      } else if (dist == "studmult") {
        .permtraces.Studmult.Hier(x, dev, col)
      }
    }
    if (lik %in% c(1, 2)) {
      ## log ##
      .permtraces.Log(x, dev, col)
    }
  }
)

#' Plot histograms of the parameters and weights
#' 
#' @description 
#' Calling [plotHist()] plots histograms of the sampled parameters and weights 
#' from MCMC sampling.More specifically, all component parameters, `K-1` of the 
#' weights and the posterior parameters are considered in the histogram plots. 
#' 
#' Note, this method is so far only implemented for Poisson and Binomial 
#' mixture distributions.
#' 
#' @param x An `mcmcoutputpermfixhier` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Histograms of the MCMC samples.
#' @exportMethod plotHist
#' 
#' @examples 
#' \dontrun{
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2, 
#'                  indicfix = TRUE)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc(storepost = FALSE)
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' f_outputperm <- mcmcpermute(f_output)
#' plotHist(f_outputperm)
#' }
#' 
#' @seealso 
#' * [mixturemcmc()] for performing MCMC sampling
#' * [mcmcpermute()] for permuting MCMC samples
#' * [plotTraces()] for plotting the traces of sampled values
#' * [plotDens()] for plotting densities of sampled values
#' * [plotSampRep()] for plotting sampling representations of sampled values
#' * [plotPointProc()] for plotting point processes for sampled values
#' * [plotPostDens()] for plotting the posterior density of component parameters
setMethod(
  "plotHist", signature(
    x = "mcmcoutputpermfixhier",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    dist <- x@model@dist
    if (dist == "poisson") {
      .permhist.Poisson.Hier(x, dev)
    } else if (dist == "binomial") {
      .permhist.Binomial(x, dev)
    }
  }
)

#' Plot densities of the parameters and weights
#' 
#' @description 
#' Calling [plotDens()] plots densities of the sampled parameters and weights 
#' from MCMC sampling.More specifically, all component parameters, `K-1` of the 
#' weights and the posterior parameters are considered in the density plots. 
#' 
#' Note, this method is so far only implemented for mixtures of Poisson or 
#' Binomial distributions.
#' 
#' @param x An `mcmcoutputpermfixhier` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Densities of the MCMC samples.
#' @exportMethod plotDens
#' 
#' @examples 
#' \dontrun{
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2, 
#'                  indicfix = TRUE)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc(storepost = FALSE)
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' f_outputperm <- mcmcpermute(f_output)
#' plotDens(f_outputperm)
#' }
#' 
#' @seealso 
#' * [mixturemcmc()] for performing MCMC sampling
#' * [mcmcpermute()] for permuting MCMC samples
#' * [plotTraces()] for plotting the traces of sampled values
#' * [plotHist()] for plotting histograms of sampled values
#' * [plotSampRep()] for plotting sampling representations of sampled values
#' * [plotPointProc()] for plotting point processes for sampled values
#' * [plotPostDens()] for plotting the posterior density of component parameters
setMethod(
  "plotDens", signature(
    x = "mcmcoutputpermfixhier",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    dist <- x@model@dist
    if (dist == "poisson") {
      .permdens.Poisson.Hier(x, dev)
    } else if (dist == "binomial") {
      .permdens.Binomial(x, dev)
    }
  }
)

#' Plot point processes of the component parameters
#' 
#' @description 
#' Calling [plotPointProc()] plots point processes of the sampled component 
#' parameters from MCMC sampling.  
#' 
#' Note, this method is so far only implemented for mixture models of Poisson 
#' or Binomial distributons.
#' 
#' @param x An `mcmcoutputpermfixhier` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Point process of the MCMC samples.
#' @exportMethod plotPointProc
#' 
#' @examples 
#' \dontrun{
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2, 
#'                  indicfix = TRUE)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc(storepost = FALSE)
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' f_outputperm <- mcmcpermute(f_output)
#' plotPointProc(f_outputperm)
#' }
#' 
#' @seealso 
#' * [mixturemcmc()] for performing MCMC sampling
#' * [mcmcpermute()] for permuting MCMC samples
#' * [plotTraces()] for plotting the traces of sampled values
#' * [plotHist()] for plotting histograms of sampled values
#' * [plotDens()] for plotting densities of sampled values
#' * [plotSampRep()] for plotting sampling representations of sampled values
#' * [plotPostDens()] for plotting posterior densities for sampled values
setMethod(
  "plotPointProc", signature(
    x = "mcmcoutputpermfixhier",
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

#' Plot sampling representations for the component parameters
#' 
#' @description 
#' Calling [plotSampRep()] plots sampling representations of the sampled 
#' component parameters from MCMC sampling.  
#' 
#' @param x An `mcmcoutputpermfixhier` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Sampling representation of the MCMC samples.
#' @exportMethod plotSampRep
#' 
#' @examples 
#' \dontrun{
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2, 
#'                  indicfix = TRUE)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc(storepost = FALSE)
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' f_outputperm <- mcmcpermute(f_output)
#' plotSampRep(f_outputperm)
#' }
#' 
#' @seealso 
#' * [mixturemcmc()] for performing MCMC sampling
#' * [mcmcpermute()] for permuting MCMC samples
#' * [plotTraces()] for plotting the traces of sampled values
#' * [plotHist()] for plotting histograms of sampled values
#' * [plotDens()] for plotting densities of sampled values
#' * [plotPointProc()] for plotting point processes of sampled values
#' * [plotPostDens()] for plotting posterior densities for sampled values
setMethod(
  "plotSampRep", signature(
    x = "mcmcoutputpermfixhier",
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

#' Plot posterior densities of the component parameters
#' 
#' @description 
#' Calling [plotPostDens()] plots posterior densities of the sampled component 
#' parameters from MCMC sampling, if the number of components is two. 
#' 
#' Note, this method is so far only implemented for Poisson and Binomial 
#' mixture distributions.
#' 
#' @param x An `mcmcoutputpermfixhier` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Posterior densities of the MCMC samples.
#' @exportMethod plotPostDens
#' @noRd
#' 
#' @examples 
#' \dontrun{
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2, 
#'                  indicfix = TRUE)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc(storepost = FALSE)
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' f_outputperm <- mcmcpermute(f_output)
#' plotPostDens(f_outputperm)
#' }
#' 
#' @seealso 
#' * [mixturemcmc()] for performing MCMC sampling
#' * [mcmcpermute()] for permuting MCMC samples
#' * [plotTraces()] for plotting the traces of sampled values
#' * [plotHist()] for plotting histograms of sampled values
#' * [plotDens()] for plotting densities of sampled values
#' * [plotSampRep()] for plotting sampling representations of sampled values
#' * [plotPointProc()] for plotting point processes for sampled values
setMethod(
  "plotPostDens", signature(
    x = "mcmcoutputpermfixhier",
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

### Traces

#' Plots traces of Poisson mixture samples
#' 
#' @description 
#' For internal usage only. This function plots the traces for sampled values 
#' from a Poisson mixture model, if an hierarchical prior has been used in 
#' sampling. The hyperparameter `b` of the hierarchical Gamma distribution is
#' plotted next to the component parameter traces.
#' 
#' @param x An `mcmcoutputpermfixhier` object containing all samples.
#' @param dev A logical indicating if the plot should be shown by a graphical 
#'   device.
#' @return A plot of the traces of sampled values.
#' @noRd
#' 
#' @seealso 
#' * [plotTraces()] for the calling function
".permtraces.Poisson.Hier" <- function(x, dev) {
  K <- x@model@K
  trace.n <- K + 1
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
    axis(2, las = 2, cex.axis = 0.7)
    mtext(
      side = 2, las = 2, bquote(lambda[k = .(k)]),
      cex = 0.6, line = 3
    )
  }
  b <- x@hyperperm$b
  plot(b,
    type = "l", axes = F,
    col = "gray68", xlab = "", ylab = ""
  )
  axis(2, las = 2, cex.axis = 0.7)
  mtext(side = 2, las = 2, "b", cex = 0.6, line = 3)
  axis(1)
  mtext(side = 1, "Iterations", cex = 0.7, line = 3)
}

#' Plots traces of Normal mixture samples
#' 
#' @description 
#' For internal usage only. This function plots the traces for sampled values 
#' from a Normal mixture model. The parameters of the hierarchical prior are 
#' plotted together with the component parameters.
#' 
#' @param x An `mcmcoutputpermfixhier` object containing all samples.
#' @param dev A logical indicating if the plot should be shown by a grapical 
#'   device.
#' @return A plot of the traces of sampled values.
#' @noRd
#' 
#' @seealso 
#' * [plotTraces()] for the calling function
".permtraces.Normal.Hier" <- function(x, dev) {
  K <- x@model@K
  trace.n <- 2 * K + 1
  if (.check.grDevice() && dev) {
    dev.new(title = "Traceplots")
  }
  par(
    mfrow = c(trace.n, 1), mar = c(1, 0, 0, 0),
    oma = c(4, 5, 4, 4)
  )
  mu <- x@parperm$mu
  sigma <- x@parperm$sigma
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
  C <- x@hyperperm$C
  plot(c,
    type = "l", axes = F,
    col = "gray68", xlab = "", ylab = ""
  )
  axis(2, las = 2, cex.axis = 0.7)
  mtext(side = 2, las = 2, "C", cex = .6, line = 3)
  axis(1)
  mtext(side = 1, "Iterations", cex = .7, line = 3)
}

#' Plots traces of Student-t mixture samples
#' 
#' @description 
#' For internal usage only. This function plots the traces for sampled values 
#' from a Student-t mixture model. The parameters of the hierarchical prior 
#' are plotted together with the component parameters.
#' 
#' @param x An `mcmcoutputpermfixhier` object containing all samples.
#' @param dev A logical indicating if the plot should be shown by a graphical 
#'   device.
#' @return A plot of the traces of sampled values.
#' @noRd
#' 
#' @seealso 
#' * [plotTraces()] for the calling function
".permtraces.Student.Hier" <- function(x, dev) {
  K <- x@model@K
  trace.n <- 3 * K + 1
  if (.check.grDevice() && dev) {
    dev.new(title = "Traceplots")
  }
  par(
    mfrow = c(trace.n, 1), mar = c(1, 0, 0, 0),
    oma = c(4, 5, 4, 4)
  )
  mu <- x@parperm$mu
  sigma <- x@parperm$sigma
  df <- x@parperm$df
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
  C <- x@hyperperm$C
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

#' Plots traces of multivariate normal mixture samples
#' 
#' @description 
#' For internal usage only. This function plots the traces for sampled values 
#' from a multivariate normal mixture model. The parameters of the hierarchical 
#' prior are plotted alongside the component parameters.
#' 
#' @param x An `mcmcoutputpermfixhier` object containing all samples.
#' @param dev A logical indicating if the plot should be shown by a grapical 
#'   device.
#' @return A plot of the traces of sampled values.
#' @noRd
#' 
#' @seealso 
#' * [plotTraces()] for the calling function
".permtraces.Normult.Hier" <- function(x, dev, col) {
  .permtraces.Normult(x, dev, col)
  r <- x@model@r
  K <- x@model@K
  C <- x@hyperperm$C
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

#' Plots traces of multivariate Student-t mixture samples
#' 
#' @description 
#' For internal usage only. This function plots the traces for sampled values 
#' from a multivariate Student-t mixture model. The parameters of the 
#' hierarchical prior are plotted alongside the component parameters.
#' 
#' @param x An `mcmcoutputpermfixhier` object containing all samples.
#' @param dev A logical indicating if the plot should be shown by a graphical 
#'   device.
#' @return A plot of the traces of sampled values.
#' @noRd
#' 
#' @seealso 
#' * [plotTraces()] for the calling function
".permtraces.Studmult.Hier" <- function(x, dev, col) {
  .permtraces.Studmult(x, dev, col)
  r <- x@model@r
  K <- x@model@K
  C <- x@hyperperm$C
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

### Histograms
### Histograms Poisson: Plots histograms for all Poisson
### parameters and the hyper-parameter 'b'.
#' Plot histograms of Poisson samples 
#' 
#' @description 
#' For internal usage only. This function plots histograms of sampled Poisson 
#' parameters and weights. In addition the parameters of the hierarchical prior 
#' `b` are plotted.
#' 
#' @param x An `mcmcoutputpermfixhier` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with histograms for the sampled parameters and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotHist()] for the calling function
".permhist.Poisson.Hier" <- function(x, dev) {
  K <- x@model@K
  if (.check.grDevice() && dev) {
    dev.new(title = "Histograms (permuted)")
  }
  lambda <- x@parperm$lambda
  b <- x@hyperperm$b
  vars <- cbind(lambda, b)
  lab.names <- vector("list", K + 1)
  for (k in 1:K) {
    lab.names[[k]] <- bquote(lambda[.(k)])
  }
  lab.names[[K + 1]] <- "b"
  .symmetric.Hist(vars, lab.names)
}

### Densities

#' Plot densities of Poisson samples
#' 
#' @description 
#' For internal usage only. This function plots densities of sampled Poisson 
#' parameters and weights. In addition the parameters of the hierarchical prior 
#' are plotted.
#' 
#' @param x An `mcmcoutputpermfixhier` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with densities for the sampled parameters and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotDens()] for the calling function
".permdens.Poisson.Hier" <- function(x, dev) {
  K <- x@model@K
  if (.check.grDevice() && dev) {
    dev.new(title = "Histograms (permuted)")
  }
  lambda <- x@parperm$lambda
  b <- x@hyperperm$b
  vars <- cbind(lambda, b)
  lab.names <- vector("list", K + 1)
  for (k in 1:K) {
    lab.names[[k]] <- bquote(lambda[.(k)])
  }
  lab.names[[K + 1]] <- "b"
  .symmetric.Dens(vars, lab.names)
}