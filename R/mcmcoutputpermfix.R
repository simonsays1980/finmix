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

#' Finmix `mcmcoutputpermfix` class 
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
#' Note that this class inherits all of its slots from the parent classes.
#' 
#' @exportClass mcmcoutputpermfix
#' @rdname mcmcoutputpermfix-class
#' @seealso 
#' * [mcmcoutputfix-class] for the parent class
#' * [mcmcpermfix-class] for the parent class
#' * [mcmcpermute()] for performing permutation of MCMC samples
.mcmcoutputpermfix <- setClass("mcmcoutputpermfix",
  contains = c("mcmcpermfix", "mcmcoutputfix"),
  validity = function(object) {
    ## else: OK
    TRUE
  }
)

#' Initializer of the `mcmcoutputpermfix` class
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
#' 
#' @keywords internal
#' 
#' @seealso 
#' * [Classes_Details] for details of class definitions, and 
#' * [setOldClass] for the relation to S3 classes
setMethod(
  "initialize", "mcmcoutputpermfix",
  function(.Object, mcmcoutput, Mperm = integer(),
           parperm = list(), logperm = list()) {
    .Object@M <- mcmcoutput@M
    .Object@burnin <- mcmcoutput@burnin
    .Object@ranperm <- mcmcoutput@ranperm
    .Object@par <- mcmcoutput@par
    .Object@log <- mcmcoutput@log
    .Object@model <- mcmcoutput@model
    .Object@prior <- mcmcoutput@prior
    .Object@Mperm <- Mperm
    .Object@parperm <- parperm
    .Object@logperm <- logperm
    .Object
  }
)

#' Shows a summary of an `mcmcoutputpermfix` object.
#' 
#' @description
#' Calling [show()] on an `mcmcoutputpermfix` object gives an overview 
#' of the `mcmcoutputpermfix` object.
#' 
#' @param object An `mcmcoutputpermfix` object.
#' @returns A console output listing the slots and summary information about
#'   each of them. 
#' @exportMethod show
#' @noRd
setMethod(
  "show", "mcmcoutputpermfix",
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
#' @param x An `mcmcoutputpermfix` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param lik An integer indicating, if the log-likelihood traces should be 
#'   plotted (default). If set to `0` the traces for the parameters 
#'   and weights are plotted instead. 
#' @param col A logical indicating, if the plot should be colored.
#' @param ... Further arguments to be passed to the plotting function.
#' @return A plot of the traces of the MCMC samples.
#' @exportMethod plotTraces 
#' @describeIn mcmcoutputpermfix-class
#' 
#' @examples 
#' # Define a Poisson mixture model with two components.
#' f_model <- model("poisson", par = list(lambda = c(0.3, 1.2)), K = 2, 
#'                  indicfix = TRUE)
#' # Simulate data from the mixture model.
#' f_data <- simulate(f_model)
#' # Define the hyper-parameters for MCMC sampling.
#' f_mcmc <- mcmc(storepost = FALSE)
#' # Define the prior distribution by relying on the data.
#' f_prior <- priordefine(f_data, f_model)
#' # Do not use a hierarchical prior.
#' setHier(f_prior) <- FALSE
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' f_outputperm <- mcmcpermute(f_output)
#' plotTraces(f_outputperm, lik = 0)
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
    x = "mcmcoutputpermfix",
    dev = "ANY",
    lik = "ANY",
    col = "ANY"
  ),
  function(x, dev = TRUE, lik = 1, col = FALSE, ...) {
    dist <- x@model@dist
    if (lik %in% c(0, 1)) {
      if (dist == "poisson") {
        .permtraces.Poisson(x, dev)
      } else if (dist == "binomial") {
        .permtraces.Binomial(x, dev)
      } else if (dist == "exponential") {
        .permtraces.Exponential(x, dev)
      } else if (dist == "normal") {
        .permtraces.Normal(x, dev)
      } else if (dist == "student") {
        .permtraces.Student(x, dev)
      } else if (dist == "normult") {
        .permtraces.Normult(x, dev, col)
      } else if (dist == "studmult") {
        .permtraces.Studmult(x, dev, col)
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
#' Calling [plotHist()] plots histograms of the sampled component parameters  
#' from MCMC sampling. 
#' 
#' @param x An `mcmcoutputpermfix` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Histograms of the MCMC samples.
#' @exportMethod plotHist 
#' @describeIn mcmcoutputpermfix-class
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
#' # Do not use a hierarchical prior.
#' setHier(f_prior) <- FALSE
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
    x = "mcmcoutputpermfix",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    dist <- x@model@dist
    if (dist == "poisson") {
      .permhist.Poisson(x, dev)
    } else if (dist == "binomial") {
      .permhist.Binomial(x, dev)
    }
  }
)

#' Plot densities of the parameters and weights
#' 
#' @description 
#' Calling [plotDens()] plots densities of the sampled component parameters  
#' from MCMC sampling. 
#' 
#' @param x An `mcmcoutputpermfix` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Densities of the MCMC samples.
#' @exportMethod plotDens
#' @describeIn mcmcoutputpermfix-class
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
#' # Do not use a hierarchical prior.
#' setHier(f_prior) <- FALSE
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
    x = "mcmcoutputpermfix",
    dev = "ANY"
  ),
  function(x, dev = TRUE, ...) {
    dist <- x@model@dist
    if (dist == "poisson") {
      .permdens.Poisson(x, dev)
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
#' @param x An `mcmcoutputpermfix` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Densities of the MCMC samples.
#' @exportMethod plotPointProc
#' @describeIn mcmcoutputpermfix-class
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
#' # Do not use a hierarchical prior.
#' setHier(f_prior) <- FALSE
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
    x = "mcmcoutputpermfix",
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
#' @param x An `mcmcoutputpermfix` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Sampling represetnation of the MCMC samples.
#' @exportMethod plotSampRep
#' @describeIn mcmcoutputpermfix-class
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
#' # Do not use a hierarchical prior.
#' setHier(f_prior) <- FALSE
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
    x = "mcmcoutputpermfix",
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
#' @param x An `mcmcoutputpermfix` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Posterior densities of the MCMC samples.
#' @exportMethod plotPostDens
#' @describeIn mcmcoutputpermfix-class
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
#' # Do not use a hierarchical prior.
#' setHier(f_prior) <- FALSE
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
    x = "mcmcoutputpermfix",
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
#' from a Poisson mixture model.
#' 
#' @param x An `mcmcoutputpermfix` object containing all samples.
#' @param dev A logical indicating if the plot should be shown by a graphical 
#'   device.
#' @return A plot of the traces of sampled values.
#' @noRd
#' 
#' @seealso 
#' * [plotTraces()] for the calling function
".permtraces.Poisson" <- function(x, dev) {
  K <- x@model@K
  trace.n <- K
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
  axis(1)
  mtext(side = 1, "Iterations", cex = 0.7, line = 3)
}

#' Plots traces of Binomial mixture samples
#' 
#' @description 
#' For internal usage only. This function plots the traces for sampled values 
#' from a Binomial mixture model.
#' 
#' @param x An `mcmcoutputpermfix` object containing all samples.
#' @param dev A logical indicating if the plot should be shown by a grapical 
#'   device.
#' @return A plot of the traces of sampled values.
#' @noRd
#' 
#' @seealso 
#' * [plotTraces()] for the calling function
".permtraces.Binomial" <- function(x, dev) {
  K <- x@model@K
  trace.n <- K
  if (.check.grDevice() && dev) {
    dev.new(title = "Traceplots")
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
  axis(1)
  mtext(side = 1, "Iterations", cex = 0.7, line = 3)
}

#' Plots traces of exponential mixture samples
#' 
#' @description 
#' For internal usage only. This function plots the traces for sampled values 
#' from a exponential mixture model.
#' 
#' @param x An `mcmcoutputpermfix` object containing all samples.
#' @param dev A logical indicating if the plot should be shown by a grapical 
#'   device.
#' @return A plot of the traces of sampled values.
#' @noRd
#' 
#' @seealso 
#' * [plotTraces()] for the calling function
".permtraces.Exponential" <- function(x, dev) {
  K <- x@model@K
  trace.n <- K
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
  axis(1)
  mtext(side = 1, "Iterations", cex = .7, line = 3)
}

#' Plots traces of Normal mixture samples
#' 
#' @description 
#' For internal usage only. This function plots the traces for sampled values 
#' from a Normal mixture model.
#' 
#' @param x An `mcmcoutputpermfix` object containing all samples.
#' @param dev A logical indicating if the plot should be shown by a grapical 
#'   device.
#' @return A plot of the traces of sampled values.
#' @noRd
#' 
#' @seealso 
#' * [plotTraces()] for the calling function
".permtraces.Normal" <- function(x, dev) {
  K <- x@model@K
  trace.n <- 2 * K
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
  axis(1)
  mtext(side = 1, "Iterations", cex = .7, line = 3)
}

#' Plots traces of Student-t mixture samples
#' 
#' @description 
#' For internal usage only. This function plots the traces for sampled values 
#' from a Student-t mixture model.
#' 
#' @param x An `mcmcoutputpermfix` object containing all samples.
#' @param dev A logical indicating if the plot should be shown by a grapical 
#'   device.
#' @return A plot of the traces of sampled values.
#' @noRd
#' 
#' @seealso 
#' * [plotTraces()] for the calling function
".permtraces.Student" <- function(x, dev) {
  K <- x@model@K
  trace.n <- 3 * K
  if (.check.grDevice() && dev) {
    dev.new(title = "Traceplots")
  }
  par(
    mfrow = c(trace.n, 1), mar = c(1, 0, 0, 0),
    oma = c(4, 5, 4, 4)
  )
  mu <- x@parperm$mu
  sigma <- x@parperm$sigma
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
  axis(1)
  mtext(side = 1, "Iterations", cex = .7, line = 3)
}

#' Plots traces of multivariate normal mixture samples
#' 
#' @description 
#' For internal usage only. This function plots the traces for sampled values 
#' from a multivariate normal mixture model.
#' 
#' @param x An `mcmcoutputpermfix` object containing all samples.
#' @param dev A logical indicating if the plot should be shown by a grapical 
#'   device.
#' @return A plot of the traces of sampled values.
#' @noRd
#' 
#' @seealso 
#' * [plotTraces()] for the calling function
".permtraces.Normult" <- function(x, dev, col) {
  K <- x@model@K
  r <- x@model@r
  if (.check.grDevice() && dev) {
    dev.new(title = "Traceplots")
  }
  trace.n <- r + 2
  par(
    mfrow = c(trace.n, 1), mar = c(1, 2, 0, 0),
    oma = c(4, 5, 2, 4)
  )
  mu <- x@parperm$mu
  sigma <- x@parperm$sigma
  if (col) {
    cscale <- rainbow(K, start = 0.5, end = 0)
  } else {
    cscale <- gray.colors(K, start = 0.5, end = 0.15)
  }
  for (rr in 1:r) {
    mmax <- max(mu[, rr, ])
    mmin <- min(mu[, rr, ])
    plot(mu[, rr, 1],
      type = "l", axes = F,
      col = cscale[1], xlab = "", ylab = "",
      ylim = c(mmin, mmax + 0.3 * (mmax - mmin))
    )
    for (k in 2:K) {
      lines(mu[, rr, k], col = cscale[k])
    }
    axis(2, las = 2, cex.axis = .7)
    mtext(
      side = 2, las = 2, bquote(mu[rr = .(rr)]),
      cex = .6, line = 3
    )
    if (rr == 1) {
      name <- vector("character", K)
      for (k in 1:K) {
        name[k] <- paste("k = ", k, sep = "")
      }
      legend("top",
        legend = name, col = cscale, horiz = TRUE,
        lty = 1
      )
    }
  }
  sigma.tr <- array(numeric(), dim = c(x@M, K))
  sigma.det <- array(numeric(), dim = c(x@M, K))
  for (k in 1:K) {
    sigma.tr[, k] <- sapply(
      seq(1, x@M),
      function(i) sum(diag(qinmatr(sigma[i, , k])))
    )
    sigma.det[, k] <- sapply(
      seq(1, x@M),
      function(i) log(det(qinmatr(sigma[i, , k])))
    )
  }
  # Sigma traces
  mmax <- max(sigma.tr)
  mmin <- min(sigma.tr)
  plot(sigma.tr[, 1],
    type = "l", axes = F,
    col = cscale[1], xlab = "", ylab = "",
    ylim = c(mmin, mmax)
  )
  for (k in 2:K) {
    lines(sigma.tr[, k], col = cscale[k])
  }
  axis(2, las = 2, cex.axis = .7)
  mtext(
    side = 2, las = 2, bquote(tr(Sigma)),
    cex = .6, line = 3
  )

  # Sigma logdets
  mmax <- max(sigma.det)
  mmin <- min(sigma.det)
  plot(sigma.det[, 1],
    type = "l", axes = F,
    col = cscale[1], xlab = "", ylab = "",
    ylim = c(mmin, mmax)
  )
  for (k in 2:K) {
    lines(sigma.det[, k], col = cscale[k])
  }
  axis(2, las = 2, cex.axis = .7)
  mtext(
    side = 2, las = 2, bquote(log(det(Sigma))),
    cex = .6, line = 3
  )

  axis(1)
  mtext(side = 1, "Iterations", cex = .7, line = 3)

  # Get moments
  moms <- permmoments_cc(x)
  for (rr in 1:r) {
    if (.check.grDevice() && dev) {
      dev.new(title = paste("Traceplots Feature ", rr, sep = ""))
    }
    par(
      mfrow = c(2, 2), mar = c(4, 4, 0.5, 0.5),
      oma = c(1.5, 2, 1, 1)
    )
    # Mu
    plot(moms$mean[, rr],
      type = "l", axes = F,
      xlab = "", ylab = "", col = cscale[K]
    )
    axis(2, las = 2, cex.axis = .7)
    mtext(
      side = 2, las = 2, bquote(mu), cex = .6,
      line = 3
    )
    axis(1)
    mtext(side = 1, "Iterations", cex = .7, line = 3)
    # Variance
    plot(moms$var[, rr],
      type = "l", axes = F,
      xlab = "", ylab = "", col = cscale[K]
    )
    axis(2, las = 2, cex.axis = .7)
    mtext(
      side = 2, las = 2, bquote(sigma), cex = .6,
      line = 3
    )
    axis(1)
    mtext(side = 1, "Iterations", cex = .7, line = 3)
    # Skewness
    plot(moms$skewness[, rr],
      type = "l", axes = F,
      xlab = "", ylab = "", col = cscale[K]
    )
    axis(2, las = 2, cex.axis = .7)
    mtext(
      side = 2, "Skewness", cex = .6,
      line = 3
    )
    axis(1)
    mtext(side = 1, "Iterations", cex = .7, line = 3)
    # Kurtosis
    plot(moms$kurtosis[, rr],
      type = "l", axes = F,
      xlab = "", ylab = "", col = cscale[K]
    )
    axis(2, las = 2, cex.axis = .7)
    mtext(
      side = 2, "Kurtosis", cex = .6,
      line = 3
    )
    axis(1)
    mtext(side = 1, "Iterations", cex = .7, line = 3)
  }
}

#' Plots traces of multivariate Student-t mixture samples
#' 
#' @description 
#' For internal usage only. This function plots the traces for sampled values 
#' from a multivariate Student-t mixture model.
#' 
#' @param x An `mcmcoutputpermfix` object containing all samples.
#' @param dev A logical indicating if the plot should be shown by a grapical 
#'   device.
#' @return A plot of the traces of sampled values.
#' @noRd
#' 
#' @seealso 
#' * [plotTraces()] for the calling function
".permtraces.Studmult" <- function(x, dev, col) {
  K <- x@model@K
  r <- x@model@r
  if (.check.grDevice() && dev) {
    dev.new(title = "Traceplots")
  }
  trace.n <- r + 2
  par(
    mfrow = c(trace.n, 1), mar = c(1, 2, 0, 0),
    oma = c(4, 5, 4, 4)
  )
  mu <- x@parperm$mu
  sigma <- x@parperm$sigma
  if (col) {
    cscale <- rainbow(K, start = 0.5, end = 0)
  } else {
    cscale <- gray.colors(K, start = 0.5, end = 0.15)
  }
  for (rr in 1:r) {
    mmax <- max(mu[, rr, ])
    mmin <- min(mu[, rr, ])
    plot(mu[, rr, 1],
      type = "l", axes = F,
      col = cscale[1], xlab = "", ylab = "",
      ylim = c(mmin, mmax + 0.3 * (mmax - mmin))
    )
    for (k in 2:K) {
      lines(mu[, rr, k], col = cscale[k])
    }
    axis(2, las = 2, cex.axis = .7)
    mtext(
      side = 2, las = 2, bquote(mu[rr = .(rr)]),
      cex = .6, line = 3
    )
    if (rr == 1) {
      name <- vector("character", K)
      for (k in 1:K) {
        name[k] <- paste("k = ", k, sep = "")
      }
      legend("top",
        legend = name, col = cscale, horiz = TRUE,
        lty = 1
      )
    }
  }
  sigma.tr <- array(numeric(), dim = c(x@M, K))
  sigma.det <- array(numeric(), dim = c(x@M, K))
  for (k in 1:K) {
    sigma.tr[, k] <- sapply(
      seq(1, x@M),
      function(i) sum(diag(qinmatr(sigma[i, , k])))
    )
    sigma.det[, k] <- sapply(
      seq(1, x@M),
      function(i) log(det(qinmatr(sigma[i, , k])))
    )
  }
  # Sigma traces
  mmax <- max(sigma.tr)
  mmin <- min(sigma.tr)
  plot(sigma.tr[, 1],
    type = "l", axes = F,
    col = cscale[1], xlab = "", ylab = "",
    ylim = c(mmin, mmax)
  )
  for (k in 2:K) {
    lines(sigma.tr[, k], col = cscale[k])
  }
  axis(2, las = 2, cex.axis = .7)
  mtext(
    side = 2, las = 2, bquote(tr(Sigma)),
    cex = .6, line = 3
  )

  # Sigma logdets
  mmax <- max(sigma.det)
  mmin <- min(sigma.det)
  plot(sigma.det[, 1],
    type = "l", axes = F,
    col = cscale[1], xlab = "", ylab = "",
    ylim = c(mmin, mmax)
  )
  for (k in 2:K) {
    lines(sigma.det[, k], col = cscale[k])
  }
  axis(2, las = 2, cex.axis = .7)
  mtext(
    side = 2, las = 2, bquote(log(det(Sigma))),
    cex = .6, line = 3
  )

  axis(1)
  mtext(side = 1, "Iterations", cex = .7, line = 3)

  # Degrees of freedom
  if (.check.grDevice() && dev) {
    dev.new(title = "Traceplots Degrees of Freedom")
  }
  degf <- x@parperm$df
  par(
    mfrow = c(K, 1), mar = c(1, 2, 0, 0),
    oma = c(4, 5, 4, 4)
  )
  for (k in 1:K) {
    plot(degf[, k],
      type = "l", axes = F,
      col = cscale[K], xlab = "", ylab = ""
    )
    axis(2, las = 2, cex.axis = .7)
    mtext(
      side = 2, las = 2, bquote(nu[k = .(k)]),
      cex = .6, line = 3
    )
  }
  axis(1)
  mtext(side = 1, "Iterations", cex = .7, line = 3)
  # Get moments
  moms <- permmoments_cc(x)
  for (rr in 1:r) {
    if (.check.grDevice() && dev) {
      dev.new(title = paste("Traceplots Feature ", rr, sep = ""))
    }
    par(
      mfrow = c(2, 2), mar = c(4, 4, 0.5, 0.5),
      oma = c(1.5, 2, 1, 1)
    )
    # Mu
    plot(moms$mean[, rr],
      type = "l", axes = F,
      xlab = "", ylab = "", col = cscale[K]
    )
    axis(2, las = 2, cex.axis = .7)
    mtext(
      side = 2, las = 2, bquote(mu), cex = .6,
      line = 3
    )
    axis(1)
    mtext(side = 1, "Iterations", cex = .7, line = 3)
    # Variance
    plot(moms$var[, rr],
      type = "l", axes = F,
      xlab = "", ylab = "", col = cscale[K]
    )
    axis(2, las = 2, cex.axis = .7)
    mtext(
      side = 2, las = 2, bquote(sigma), cex = .6,
      line = 3
    )
    axis(1)
    mtext(side = 1, "Iterations", cex = .7, line = 3)
    # Skewness
    plot(moms$skewness[, rr],
      type = "l", axes = F,
      xlab = "", ylab = "", col = cscale[K]
    )
    axis(2, las = 2, cex.axis = .7)
    mtext(
      side = 2, "Skewness", cex = .6,
      line = 3
    )
    axis(1)
    mtext(side = 1, "Iterations", cex = .7, line = 3)
    # Kurtosis
    plot(moms$kurtosis[, rr],
      type = "l", axes = F,
      xlab = "", ylab = "", col = cscale[K]
    )
    axis(2, las = 2, cex.axis = .7)
    mtext(
      side = 2, "Kurtosis", cex = .6,
      line = 3
    )
    axis(1)
    mtext(side = 1, "Iterations", cex = .7, line = 3)
  }
}

#' Plots traces of log-likelihood samples
#' 
#' @description 
#' For internal usage only. This function plots the traces for sampled values 
#' from any mixture model.
#' 
#' @param x An `mcmcoutputpermfix` object containing all samples.
#' @param dev A logical indicating if the plot should be shown by a grapical 
#'   device.
#' @return A plot of the traces of sampled values.
#' @noRd
#' 
#' @seealso 
#' * [plotTraces()] for the calling function
".permtraces.Log" <- function(x, dev, col) {
  if (.check.grDevice() && dev) {
    dev.new(title = "Log Likelihood Traceplots")
  }
  if (col) {
    cscale <- rainbow(3, start = 0, end = .5)
  } else {
    cscale <- gray.colors(3, start = 0, end = .15)
  }
  par(
    mfrow = c(2, 1), mar = c(1, 0, 0, 0),
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
  axis(1)
  mtext(side = 1, "Iterations", cex = 0.7, line = 3)
}

### Histograms

#' Plot histograms of Poisson samples
#' 
#' @description 
#' For internal usage only. This function plots histograms of sampled Poisson 
#' parameters.
#' 
#' @param x An `mcmcoutputpermfix` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with histograms for the sampled parameters and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotHist()] for the calling function
".permhist.Poisson" <- function(x, dev) {
  K <- x@model@K
  if (.check.grDevice() && dev) {
    dev.new(title = "Histograms (permuted)")
  }
  lambda <- x@parperm$lambda
  lab.names <- vector("list", K)
  for (k in 1:K) {
    lab.names[[k]] <- bquote(lambda[.(k)])
  }
  .symmetric.Hist(lambda, lab.names)
}

#' Plot histograms of Binomial samples
#' 
#' @description 
#' For internal usage only. This function plots histograms of sampled Binomial 
#' parameters.
#' 
#' @param x An `mcmcoutputpermfix` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with histograms for the sampled parameters and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotHist()] for the calling function
".permhist.Binomial" <- function(x, dev) {
  K <- x@model@K
  if (.check.grDevice() && dev) {
    dev.new(title = "Histograms (permuted)")
  }
  p <- x@parperm$p
  lab.names <- vector("list", K)
  for (k in 1:K) {
    lab.names[[k]] <- bquote(p[.(k)])
  }
  .symmetric.Hist(p, lab.names)
}
### Densities

#' Plot densities of Poisson samples
#' 
#' @description 
#' For internal usage only. This function plots densities of sampled Poisson 
#' parameters.
#' 
#' @param x An `mcmcoutputpermfix` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with densities for the sampled parameters and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotDens()] for the calling function
".permdens.Poisson" <- function(x, dev) {
  K <- x@model@K
  if (.check.grDevice() && dev) {
    dev.new(title = "Densities (permuted)")
  }
  lambda <- x@parperm$lambda
  lab.names <- vector("list", K)
  for (k in 1:K) {
    lab.names[[k]] <- bquote(lambda[.(k)])
  }
  .symmetric.Dens(lambda, lab.names)
}

#' Plot densities of Binomial samples
#' 
#' @description 
#' For internal usage only. This function plots densities of sampled Binomial 
#' parameters.
#' 
#' @param x An `mcmcoutputpermfix` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with densities for the sampled parameters and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotDens()] for the calling function
".permdens.Binomial" <- function(x, dev) {
  K <- x@model@K
  if (.check.grDevice() && dev) {
    dev.new(title = "Densities (permuted)")
  }
  p <- x@parperm$p
  lab.names <- vector("list", K)
  for (k in 1:K) {
    lab.names[[k]] <- bquote(p[.(k)])
  }
  .symmetric.Dens(p, lab.names)
}

### Plot Point Processes

#' Plot point processes of Poisson samples
#' 
#' @description 
#' For internal usage only. This function plots the point process of sampled 
#' Poisson parameters and weights against a random normal sample.
#' 
#' @param x An `mcmcoutputpermfix` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with the point process for the sampled parameters and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotPointProc()] for the calling function
".permpointproc.Poisson" <- function(x, dev) {
  K <- x@model@K
  M <- x@Mperm
  if (.check.grDevice() && dev) {
    dev.new("Point Process Representation (MCMC permuted)")
  }
  # Produces an M x K grid
  y.grid <- replicate(K, rnorm(M))
  if (median(x@parperm$lambda) < 1) {
    lambda <- log(x@parperm$lambda)
  } else {
    lambda <- x@parperm$lambda
  }
  col.grid <- gray.colors(K,
    start = 0,
    end = 0.5
  )
  legend.names <- vector("list", K)
  for (k in seq(1, K)) {
    legend.names[[k]] <- bquote(lambda[.(k)])
  }
  plot(lambda, y.grid,
    pch = 20, col = col.grid,
    cex = .7, cex.axis = .7, cex.lab = .7,
    main = "", ylab = "", xlab = ""
  )
  mtext(
    side = 1, bquote(lambda), cex = .7,
    cex.lab = .7, line = 3
  )
  legend("topright",
    legend = do.call(
      expression,
      legend.names
    ),
    col = col.grid, fill = col.grid
  )
}

#' Plot point processes of Binomial samples
#' 
#' @description 
#' For internal usage only. This function plots the point process of sampled 
#' Binomial parameters and weights against a random normal sample.
#' 
#' @param x An `mcmcoutputpermfix` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with the point process for the sampled parameters and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotPointProc()] for the calling method
".permpointproc.Binomial" <- function(x, dev) {
  K <- x@model@K
  M <- x@M
  if (.check.grDevice() && dev) {
    dev.new(title = "Point Process Representation (MCMC permuted)")
  }
  y.grid <- replicate(K, rnorm(M))
  p <- x@par$p
  col.grid <- gray.colors(K,
    start = 0,
    end = 0.5
  )
  legend.names <- vector("list", K)
  for (k in seq(1, K)) {
    legend.names[[k]] <- bquote(p[.(k)])
  }
  plot(p, y.grid,
    pch = 20, col = col.grid,
    cex = .7, cex.axis = .7, cex.lab = .7,
    main = "", ylab = "", xlab = ""
  )
  mtext(
    side = 1, bquote(p), cex = .7,
    cex.lab = .7, line = 3
  )
  legend("topright",
    legend = do.call(
      expression,
      legend.names
    ),
    col = col.grid, fill = col.grid
  )
}

### Plot sampling representation

#' Plot sampling representation of Poisson samples
#' 
#' @description 
#' For internal usage only. This function plots the sampling representation of 
#' sampled Poisson parameters and weights. Each parameter sample is plotted 
#' against its permuted counterpart.
#' 
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with the sampling representation for the sampled parameters 
#'   and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotSampRep()] for the calling function
".permsamprep.Poisson" <- function(x, dev) {
  K <- x@model@K
  if (K == 1) {
    warning(paste("Sampling representation is only ",
      "available for mixture models with ",
      "K > 1.",
      sep = ""
    ))
    return(FALSE)
  }
  M <- x@Mperm
  n <- min(2000, x@Mperm)
  n.perm <- choose(K, 2) * factorial(2)
  lambda <- x@parperm$lambda
  if (.check.grDevice() && dev) {
    dev.new(title = "Sampling Representation")
  }
  comb <- as.matrix(expand.grid(seq(1, K), seq(1, K)))
  comb <- comb[which(comb[, 1] != comb[, 2]), ]
  lambda <- lambda[seq(1, n), ]
  lambda <- matrix(lambda[, comb], nrow = n * n.perm, ncol = 2)
  plot(lambda,
    col = "gray47", cex.lab = .7, cex.axis = .7,
    cex = .7, pch = 20, main = "", xlab = "", ylab = ""
  )
  abline(0, 1, lty = 1)
  mtext(
    side = 1, bquote(lambda), cex = .7, cex.lab = .7,
    line = 3
  )
  mtext(
    side = 2, bquote(lambda), cex = .7, cex.lab = .7,
    line = 3
  )
}

#' Plot sampling representation of Binomial samples
#' 
#' @description 
#' For internal usage only. This function plots the sampling representation of 
#' sampled Binomial parameters and weights. Each parameter sample is plotted 
#' against its permuted counterpart.
#' 
#' @param x An `mcmcoutputpermfix` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with the sampling representation for the sampled parameters 
#'   and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotSampRep()] for the calling function
".permsamprep.Binomial" <- function(x, dev) {
  K <- x@model@K
  if (K == 1) {
    warning(paste("Sampling representation is only ",
      "available for mixture models with ",
      "K > 1.",
      sep = ""
    ))
    return(FALSE)
  }
  M <- x@M
  n <- min(2000, x@M)
  n.perm <- choose(K, 2) * factorial(2)
  p <- x@parperm$p
  if (.check.grDevice() && dev) {
    dev.new(title = "Sampling Representation (MCMC permuted)")
  }
  comb <- as.matrix(expand.grid(seq(1, K), seq(1, K)))
  comb <- comb[which(comb[, 1] != comb[, 2]), ]
  p <- p[seq(1, n), ]
  p <- matrix(p[, comb], nrow = n * n.perm, ncol = 2)
  plot(p,
    col = "gray47", cex.lab = .7, cex.axis = .7,
    cex = .7, pch = 20, main = "", xlab = "", ylab = ""
  )
  abline(0, 1, lty = 1)
  mtext(
    side = 1, bquote(p), cex = .7, cex.lab = .7,
    line = 3
  )
  mtext(
    side = 2, bquote(p), cex = .7, cex.lab = .7,
    line = 3
  )
}

### Posterior Density

#' Plot posterior density of Poisson samples
#' 
#' @description 
#' For internal usage only. This function plots the posterior density of 
#' sampled Poisson parameters and weights. 
#' 
#' @param x An `mcmcoutputpermfix` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with the posterior density for the sampled parameters 
#'   and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotPostdens()] for the calling function
".permpostdens.Poisson" <- function(x, dev) {
  K <- x@model@K
  if (K != 2) {
    warning(paste("A plot of the posterior density is ",
      "available only for K = 2.",
      sep = ""
    ))
  } else {
    M <- x@M
    n <- min(2000, M)
    lambda <- x@parperm$lambda
    lambda <- lambda[seq(1, n), ]
    dens <- bkde2D(lambda, bandwidth = c(
      sd(lambda[, 1]),
      sd(lambda[, 2])
    ))
    if (.check.grDevice() && dev) {
      dev.new(title = "Posterior Density Contour Plot (MCMC)")
    }
    contour(dens$x1, dens$x2, dens$fhat,
      cex = .7,
      cex.lab = .7, cex.axis = .7, col = "gray47",
      main = "", xlab = "", ylab = ""
    )
    mtext(
      side = 1, bquote(lambda[1]), cex = .7,
      cex.lab = .7, line = 3
    )
    mtext(
      side = 2, bquote(lambda[2]), cex = .7,
      cex.lab = .7, line = 3
    )
    if (.check.grDevice() && dev) {
      dev.new(title = "Posterior Density Persepctive Plot (MCMC)")
    }
    persp(dens$x1, dens$x2, dens$fhat,
      col = "gray65",
      border = "gray47", theta = 55, phi = 30,
      expand = .5, lphi = 180, ltheta = 90,
      r = 40, d = .1, ticktype = "detailed", zlab =
        "Density", xlab = "k = 1", ylab = "k = 2"
    )
  }
}

#' Plot posterior density of Binomial samples
#' 
#' @description 
#' For internal usage only. This function plots the posterior density of 
#' sampled Binomial parameters and weights. 
#' 
#' @param x An `mcmcoutputpermfix` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with the posterior density for the sampled parameters 
#'   and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotPostDens()] for the calling function
".permpostdens.Binomial" <- function(x, dev) {
  K <- x@model@K
  if (K != 2) {
    warning(paste("A plot of the posterior density is ",
      "available only for K = 2.",
      sep = ""
    ))
  } else {
    M <- x@M
    n <- min(2000, M)
    p <- x@parperm$p
    p <- p[seq(1, n), ]
    dens <- bkde2D(p, bandwidth = c(
      sd(p[, 1]),
      sd(p[, 2])
    ))
    if (.check.grDevice() && dev) {
      dev.new(title = "Posterior Density Contour Plot (MCMC)")
    }
    contour(dens$x1, dens$x2, dens$fhat,
      cex = .7,
      cex.lab = .7, cex.axis = .7, col = "gray47",
      main = "", xlab = "", ylab = ""
    )
    mtext(
      side = 1, bquote(p[1]), cex = .7,
      cex.lab = .7, line = 3
    )
    mtext(
      side = 2, bquote(p[2]), cex = .7,
      cex.lab = .7, line = 3
    )
    if (.check.grDevice() && dev) {
      dev.new(title = "Posterior Density Persepctive Plot (MCMC)")
    }
    persp(dens$x1, dens$x2, dens$fhat,
      col = "gray65",
      border = "gray47", theta = 55, phi = 30,
      expand = .5, lphi = 180, ltheta = 90,
      r = 40, d = .1, ticktype = "detailed", zlab =
        "Density", xlab = "k = 1", ylab = "k = 2"
    )
  }
}