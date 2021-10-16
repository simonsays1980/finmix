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

#' Finmix `mcmcoutput` class for hierarchical priors
#' 
#' @description 
#' This class stores in addition to the information from its parent class 
#' `mcmcoutputfix` also the sampled parameters from the hierarchical prior. 
#' 
#' @slot hyper A list storing the sampled parameters from the hierarchical 
#'   prior. 
#' @exportClass mcmcoutputfixhier
#' @rdname mcmcoutputfixhier-class
#' 
#' @seealso 
#' * [mcmcoutputfix-class] for the parent class``
.mcmcoutputfixhier <- setClass("mcmcoutputfixhier",
  representation(hyper = "list"),
  contains = c("mcmcoutputfix"),
  validity = function(object) {
    ## else: OK
    TRUE
  },
  prototype(hyper = list())
)

#' Shows a summary of an `mcmcoutputfixhier` object.
#' 
#' Calling [show()] on an `mcmcoutputfixhier` object gives an overview 
#' of the `mcmcoutputfixhier` object.
#' 
#' @param object An `mcmcoutputfixhier` object.
#' @returns A console output listing the slots and summary information about
#'   each of them. 
#' @exportMethod show
#' @noRd
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

#' Plot traces of MCMC sampling
#' 
#' @description 
#' Calling [plotTraces()] plots the MCMC traces of the mixture log-likelihood 
#' , the mixture log-likelihood of the prior distribution, the log-likelihood 
#' of the complete data posterior, or the weights and parameters if `lik` is 
#' set to `0`. 
#' 
#' If `lik` is set to `0` the parameters of the components and the posterior 
#' parameters are plotted together with `K-1` weights.
#' 
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param lik An integer indicating, if the log-likelihood traces should be 
#'   plotted (default). If set to `0` the traces for the parameters 
#'   and weights are plotted instead. 
#' @param col A logical indicating, if the plot should be colored.
#' @param ... Further arguments to be passed to the plotting function.
#' @return A plot of the traces of the MCMC samples.
#' @exportMethod plotTraces
#' @noRd
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
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' plotTraces(f_output, lik = 0)
#' 
#' @seealso 
#' * [mixturemcmc()] for performing MCMC sampling
#' * [plotHist()] for plotting histograms of sampled values
#' * [plotDens()] for plotting densities of sampled values
#' * [plotSampRep()] for plotting sampling representations of sampled values
#' * [plotPointProc()] for plotting point processes for sampled values
#' * [plotPostDens()] for plotting the posterior density of component parameters
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

#' Plot histograms of the parameters and weights
#' 
#' @description 
#' Calling [plotHist()] plots histograms of the sampled parameters and weights 
#' from MCMC sampling.More specifically, all component parameters, `K-1` of the 
#' weights and the posterior parameters are considered in the histogram plots. 
#' 
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Histograms of the MCMC samples.
#' @exportMethod plotHist
#' @noRd
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
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' plotHist(f_output)
#' 
#' @seealso 
#' * [mixturemcmc()] for performing MCMC sampling
#' * [plotTraces()] for plotting the traces of sampled values
#' * [plotDens()] for plotting densities of sampled values
#' * [plotSampRep()] for plotting sampling representations of sampled values
#' * [plotPointProc()] for plotting point processes for sampled values
#' * [plotPostDens()] for plotting the posterior density of component parameters
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

#' Plot densities of the parameters and weights
#' 
#' @description 
#' Calling [plotDens()] plots densities of the sampled parameters and weights 
#' from MCMC sampling.More specifically, all component parameters, `K-1` of the 
#' weights and the posterior parameters are considered in the density plots. 
#' 
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Densities of the MCMC samples.
#' @exportMethod plotDens
#' @noRd
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
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' plotDens(f_output)
#' 
#' @seealso 
#' * [mixturemcmc()] for performing MCMC sampling
#' * [plotTraces()] for plotting the traces of sampled values
#' * [plotHist()] for plotting histograms of sampled values
#' * [plotSampRep()] for plotting sampling representations of sampled values
#' * [plotPointProc()] for plotting point processes for sampled values
#' * [plotPostDens()] for plotting the posterior density of component parameters
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

#' Plot point processes of the component parameters
#' 
#' @description 
#' Calling [plotPointProc()] plots point processes of the sampled component 
#' parameters from MCMC sampling.  
#' 
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Point process of the MCMC samples.
#' @exportMethod plotPointProc
#' @noRd
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
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' plotPointProc(f_output)
#' 
#' @seealso 
#' * [mixturemcmc()] for performing MCMC sampling
#' * [plotTraces()] for plotting the traces of sampled values
#' * [plotHist()] for plotting histograms of sampled values
#' * [plotDens()] for plotting densities of sampled values
#' * [plotSampRep()] for plotting sampling representations of sampled values
#' * [plotPostDens()] for plotting posterior densities for sampled values
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

#' Plot sampling representations for the component parameters.
#' 
#' @description 
#' Calling [plotSampRep()] plots sampling representations of the sampled 
#' component parameters from MCMC sampling.  
#' 
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Sampling representation of the MCMC samples.
#' @exportMethod plotSampRep
#' @noRd
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
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' plotSampRep(f_output)
#' 
#' @seealso 
#' * [mixturemcmc()] for performing MCMC sampling
#' * [plotTraces()] for plotting the traces of sampled values
#' * [plotHist()] for plotting histograms of sampled values
#' * [plotDens()] for plotting densities of sampled values
#' * [plotPointProc()] for plotting point processes of sampled values
#' * [plotPostDens()] for plotting posterior densities for sampled values
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

#' Plot posterior densities of the component parameters
#' 
#' @description 
#' Calling [plotPostDens()] plots posterior densities of the sampled component 
#' parameters from MCMC sampling, if the number of components is two. 
#' 
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown by a graphical 
#'   device. If plots should be stored to a file set `dev` to `FALSE`. 
#' @param ... Further arguments to be passed to the plotting function.
#' @return Posterior densities of the MCMC samples.
#' @exportMethod plotPostDens
#' @noRd
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
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' plotPostDens(f_output)
#' 
#' @seealso 
#' * [mixturemcmc()] for performing MCMC sampling
#' * [plotTraces()] for plotting the traces of sampled values
#' * [plotHist()] for plotting histograms of sampled values
#' * [plotDens()] for plotting densities of sampled values
#' * [plotSampRep()] for plotting sampling representations of sampled values
#' * [plotPointProc()] for plotting point processes for sampled values
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

#' Constructs a sub-chain of MCMC samples 
#' 
#' @description 
#' Calling [subseq()] constructs an MCMC sub-chain from the samples in the 
#' passed-in `mcmcoutput` object specfied by the index `array` in `index`. This 
#' can be advantageous, if chains are non-stationary. For successful MCMC 
#' sampling the chain must be converged to the target distribution, the true 
#' distribution of parameters, weights and indicators.
#' 
#' @param object An `mcmcoutput` object containing all sampled values.
#' @param index An array specifying the extraction of the sub-chain.
#' @return An `mcmcoutput` object containing the values from the sub-chain.
#' @exportMethod subseq
#' @noRd
setMethod(
  "subseq", signature(
    object = "mcmcoutputfixhier",
    index = "array"
  ),
  function(object, index) {
    ## Call 'subseq()' from 'mcmcoutputfix'
    object <- callNextMethod(object, index)
    dist <- object@model@dist
    ## hyper ##
    if (dist == "poisson") {
      .subseq.Poisson.Hier(object, index)
    } else if (dist %in% c("normal", "student")) {
      .subseq.Norstud.Hier(object, index)
    } else if (dist %in% c("normult", "studmult")) {
      .subseq.Normultstud.Hier(object, index)
    }
  }
)

#' Swaps elements between components
#' 
#' @description 
#' Not yet implemented.
#' 
#' @param object An `mcmcoutput` object containing the sampled values.
#' @param index An array specifying the extraction of the values.
#' @return An `mcmcoutput` object with swapped elements.
#' @exportMethod swapElements
#' @noRd 
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

#' Getter method of `mcmcoutput` class.
#' 
#' Returns the `hyper` slot.
#' 
#' @param object An `mcmcoutput` object.
#' @returns The `hyper` slot of the `object`.
#' @exportMethod getHyper
#' @noRd
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
#' # Start MCMC sampling.
#' f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
#' # Get the slot.
#' getHyper(f_output)
#' 
#' @seealso 
#' * [mcmcoutput-class] for the class definition
#' * [mixturemcmc()] for performing MCMC sampling
setMethod(
  "getHyper", "mcmcoutputfixhier",
  function(object) {
    return(object@hyper)
  }
)

## No setters for this object as it is not intended 
## that users manipulate this object. 		    	

### Private functions
### These functions are not exported.

### Plot

### Plot Traces
### Plot traces Poisson: Plots traces for each component
### parameter of a Poisson mixture and the hyper parameter 'b'.

#' Plots traces of Poisson mixture samples
#' 
#' @description 
#' For internal usage only. This function plots the traces for sampled values 
#' from a Poisson mixture model.
#' 
#' @param x An `mcmcoutput` object containing all samples.
#' @param dev A logical indicating if the plot should be shown by a grapical 
#'   device.
#' @return A plot of the traces of sampled values.
#' @noRd
#' 
#' @seealso 
#' * [plotTraces()] for the calling function
".traces.Poisson.Hier" <- function(x, dev) {
  K <- x@model@K
  trace.n <- K + 1
  if (.check.grDevice() && dev) {
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

#' Plots traces of normal mixture samples
#' 
#' @description 
#' For internal usage only. This function plots the traces for sampled values 
#' from a normal mixture model.
#' 
#' @param x An `mcmcoutput` object containing all samples.
#' @param dev A logical indicating if the plot should be shown by a grapical 
#'   device.
#' @return A plot of the traces of sampled values.
#' @noRd
#' 
#' @seealso 
#' * [plotTraces()] for the calling function
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

#' Plots traces of Student-t mixture samples
#' 
#' @description 
#' For internal usage only. This function plots the traces for sampled values 
#' from a Student-t mixture model.
#' 
#' @param x An `mcmcoutput` object containing all samples.
#' @param dev A logical indicating if the plot should be shown by a grapical 
#'   device.
#' @return A plot of the traces of sampled values.
#' @noRd
#' 
#' @seealso 
#' * [plotTraces()] for the calling function
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

#' Plots traces of multivariate normal mixture samples
#' 
#' @description 
#' For internal usage only. This function plots the traces for sampled values 
#' from a multivariate normal mixture model.
#' 
#' @param x An `mcmcoutput` object containing all samples.
#' @param dev A logical indicating if the plot should be shown by a grapical 
#'   device.
#' @return A plot of the traces of sampled values.
#' @noRd
#' 
#' @seealso 
#' * [plotTraces()] for the calling function
".traces.Normult.Hier" <- function(x, dev, col) {
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

#' Plots traces of multivariate Student-t mixture samples
#' 
#' @description 
#' For internal usage only. This function plots the traces for sampled values 
#' from a multivariate Student-t mixture model.
#' 
#' @param x An `mcmcoutput` object containing all samples.
#' @param dev A logical indicating if the plot should be shown by a grapical 
#'   device.
#' @return A plot of the traces of sampled values.
#' @noRd
#' 
#' @seealso 
#' * [plotTraces()] for the calling function
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

#' Plot histograms of Poisson samples
#' 
#' @description 
#' For internal usage only. This function plots histograms of sampled Poisson 
#' parameters and weights. In addition it plots the histogram of the 
#' parameter `b` of the hierarchical prior.
#' 
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with histograms for the smapled parameters and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotHist()] for the calling function
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

#' Plot histograms of normal samples
#' 
#' @description 
#' For internal usage only. This function plots histograms of sampled normal
#' parameters and weights. In addition it plots the sampled parameter `C` of 
#' the hierarchical prior.
#' 
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with histograms for the smapled parameters and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotHist()] for the calling function
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

#' Plot histograms of Student-t samples
#' 
#' @description 
#' For internal usage only. This function plots histograms of sampled Student-t
#' parameters and weights. In addition it plots the sampled parameter `C` of 
#' the hierarchical prior.
#' 
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with histograms for the sampled parameters and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotHist()] for the calling function
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

#' Plot histograms of multivariate normal samples
#' 
#' @description 
#' For internal usage only. This function plots histograms of sampled 
#' multivariate normal parameters and weights. In addition it plots the 
#' the logarithmised determinant and the trace of the parameter matrix `C`.
#' 
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with histograms for the smapled parameters and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotHist()] for the calling function
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

#' Plot histograms of multivariate Student-t samples
#' 
#' @description 
#' For internal usage only. This function plots histograms of sampled 
#' multivariate Student-t parameters and weights. In addition it plots the 
#' the logarithmised determinant and the trace of the parameter matrix `C`.
#' 
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with histograms for the smapled parameters and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotHist()] for the calling function
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

#' Plot densities of Poisson samples
#' 
#' @description 
#' For internal usage only. This function plots densities of sampled Poisson 
#' parameters and weights. In addition it plots the Kernel densities of the 
#' parameter `b` of the hierarchical prior.
#' 
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with densities for the sampled parameters and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotDens()] for the calling function
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

#' Plot densities of normal samples
#' 
#' @description 
#' For internal usage only. This function plots densities of sampled normal 
#' parameters and weights. In addiiton it plots the Kernel densities of the 
#' parameter `C` of the hierarchical prior.
#' 
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with densities for the sampled parameters and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotDens()] for the calling function
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

#' Plot densities of Student-t samples
#' 
#' @description 
#' For internal usage only. This function plots densities of sampled Student-t 
#' parameters and weights. In addiiton it plots the Kernel densities of the 
#' parameter `C` of the hierarchical prior.
#' 
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with densities for the sampled parameters and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotDens()] for the calling function
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

#' Plot densities of multivariate normal samples
#' 
#' @description 
#' For internal usage only. This function plots densities of sampled 
#' multivariate normal parameters and weights. In addition it plots Kernel 
#' densities of the logarithmized determinant and the trace of the parameter 
#' matrix `C` of the hierarchical prior.
#' 
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with densities for the sampled parameters and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotDens()] for the calling function
".dens.Normult.Hier" <- function(x, dev) {
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

#' Plot densities of multivariate Student-t samples
#' 
#' @description 
#' For internal usage only. This function plots densities of sampled 
#' multivariate Student-t parameters and weights. In addition it plots Kernel 
#' densities of the logarithmized determinant and the trace of the parameter 
#' matrix `C` of the hierarchical prior.
#' 
#' @param x An `mcmcoutput` object containing all sampled values.
#' @param dev A logical indicating, if the plots should be shown on a graphical 
#'   device.
#' @return A plot with densities for the sampled parameters and weights.
#' @noRd
#' 
#' @seealso 
#' * [plotDens()] for the calling function
".dens.Studmult.Hier" <- function(x, dev) {
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

#' Generates sub-chains from Poisson MCMC samples with hierarchical prior
#' 
#' @description 
#' For internal usage only. This function generates sub-chains from an 
#' `mcmcoutput` object with a Poisson `model` by defining an `index` array 
#' specifying how extraction of sub-samples should be performed. Has errors for 
#' some `mcmcoutput` sub-classes.
#' 
#' @param obj An `mcmcoutput` object containing all MCMC samples.
#' @param index An array specifying the extraction of sub-samples.
#' @return An `mcmcoutput` object containing sub-chains.
#' @noRd
#' 
#' @seealso 
#' * [subseq()] for the calling method
".subseq.Poisson.Hier" <- function(obj, index) {
  obj@hyper$b <- array(obj@hyper$b[index,],
    dim = c(obj@M, 1)
  )
  return(obj)
}

#' Generates sub-chains from Normal and Student-t MCMC samples with 
#' hierarchical prior
#' 
#' @description 
#' For internal usage only. This function generates sub-chains from an 
#' `mcmcoutput` object with a Normal or Student-t `model` by defining an 
#' `index` array specifying how extraction of sub-samples should be performed. 
#' Has errors for some `mcmcoutput` sub-classes.
#' 
#' @param obj An `mcmcoutput` object containing all MCMC samples.
#' @param index An array specifying the extraction of sub-samples.
#' @return An `mcmcoutput` object containing sub-chains.
#' @noRd
#' 
#' @seealso 
#' * [subseq()] for the calling method
".subseq.Norstud.Hier" <- function(obj, index) {
  obj@hyper$C <- array(obj@hyper$C[index,],
    dim = c(obj@M, 1)
  )
  return(obj)
}

#' Generates sub-chains from multivariate Normal and Student-t MCMC samples with 
#' hierarchical prior
#' 
#' @description 
#' For internal usage only. This function generates sub-chains from an 
#' `mcmcoutput` object with a multivariate Normal or Student-t `model` by 
#' defining an `index` array specifying how extraction of sub-samples should be 
#' performed. Has errors for some `mcmcoutput` sub-classes.
#' 
#' @param obj An `mcmcoutput` object containing all MCMC samples.
#' @param index An array specifying the extraction of sub-samples.
#' @return An `mcmcoutput` object containing sub-chains.
#' @noRd
#' 
#' @seealso 
#' * [subseq()] for the calling method
".subseq.Normultstud.Hier" <- function(obj, index) {
  obj@hyper$C <- array(obj@hyper$C[index, ],
    dim = c(obj@M, obj@model@K)
  )
  return(obj)
}