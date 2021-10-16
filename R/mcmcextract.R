#' Finmix `mcmcextract` class
#' 
#' @description
#' This is a leight-weighted class containing the major results from MCMC 
#' sampling to calculate model moments from MCMC samples. Note that momentarily
#' only methods for the multivariate Normal mixture are implemented.
#' 
#' @slot dist A character defining the finite mixture model that has been used 
#'   in MCMC sampling.
#' @slot K An integer specifying the number of components of the mixture model.
#' @slot r An integer specifying the number of dimensions of the mixture model.
#' @slot par A list storing the sample component parameters from MCMC sampling. 
#' @slot weight A n array storing the sample weight parameters from MCMC 
#' sampling.
#' 
#' @exportClass mcmcextract
#' @name mcmcextract-class
#' @noRd
.mcmcextract <- setClass("mcmcextract",
  representation(
    dist = "character",
    K = "integer",
    r = "integer",
    par = "list",
    weight = "array"
  ),
  validity = function(object) {
    ## else: OK ##
    TRUE
  },
  prototype(
    dist = character(),
    K = integer(),
    r = integer(),
    par = list(),
    weight = array()
  )
)

#' Calculate the model moments of MCMC samples
#' 
#' @description 
#' For internal usage only. This function calculates the finite mixture moments 
#' of a mixture model from the MCMC samples. Note that this function is 
#' momentarily only implemented for a mixture of multivariate Normal 
#' distributions. 
#' 
#' @param obj An `mcmcextract` object containing the parameters and weights 
#'   from MCMC sampling.
#' @return A list containing the model moments calculated from MCMC samples.
#' @exportMethod moments
#' @noRd
#' @seealso 
#' * [mcmcoutput-class] for the results from MCMC sampling
#' * [extract()] for the calling method
setMethod(
  "moments", signature(object = "mcmcextract"),
  function(object) {
    dist <- object@dist
    if (dist == "normult") {
      .moments.Normult.Mcmcextract(object)
    }
  }
)

#' Calculate the model moments of multivariate Normal MCMC samples
#' 
#' @description 
#' For internal usage only. This function calculates the finite mixture moments 
#' of a multivariate Normal mixture model from the MCMC samples. 
#' 
#' @param obj An `mcmcextract` object containing the parameters and weights 
#'   from MCMC sampling.
#' @return A list containing the model moments calculated from MCMC samples.
#' @noRd
#' @seealso 
#' * [mcmcoutput-class] for the results from MCMC sampling
#' * [extract()] for the calling method
".moments.Normult.Mcmcextract" <- function(obj) {
  K <- obj@K
  r <- obj@r
  moments.means <- apply(sapply(
    seq(length = r),
    function(i) obj@par$mu[i, ] * obj@weight
  ),
  2, sum,
  na.rm = TRUE
  )
  moments.var <- matrix(0.0, nrow = r, ncol = r)
  moments.W <- matrix(0.0, nrow = r, ncol = r)
  moments.B <- matrix(0.0, nrow = r, ncol = r)
  for (k in 1:K) {
    moments.var <- moments.var + outer(obj@par$mu[, k], obj@par$mu[, k]) +
      qinmatr(obj@par$sigma[, k]) * obj@weight[k]
    moments.W <- moments.W + qinmatr(obj@par$sigma[, k]) * obj@weight[k]
    d <- obj@par$mu[, k] - moments.means
    moments.B <- moments.B + outer(d, d) * obj@weight[k]
  }
  moments.var <- moments.var - outer(moments.means, moments.means)
  cd <- diag(1 / diag(moments.var)^0.5)
  moments.corr <- cd %*% moments.var %*% cd
  moments.Rtr <- 1 - sum(diag(moments.W)) / sum(diag(moments.var))
  moments.Rdet <- 1 - det(moments.W) / det(moments.var)
  zm <- vector("numeric", r)
  zm[seq(2, 4, by = 2)] <- exp(cumsum(log(seq(1, 4, by = 2))))
  moments.higher <- matrix(0.0, nrow = r, ncol = 4)
  for (m in 1:4) {
    for (rr in 1:r) {
      sigma <- simplify2array(sapply(seq(1, K), function(i) qinmatr(obj@par$sigma[, i]),
        simplify = FALSE
      ))[rr, rr, ]
      moments.higher[rr, m] <- sum(obj@weight * (obj@par$mu[rr, ] - moments.means[rr])^m)
      for (n in 1:m) {
        cm <- (obj@par$mu[rr, ] - moments.means[rr])^(m - n) * sigma^(n / 2) * zm[n]
        print(choose(m, n) * sum(obj@weight * cm))
        moments.higher[rr, m] <- moments.higher[rr, m] + choose(m, n) * sum(obj@weight * cm)
      }
    }
  }
  moments.skewness <- moments.higher[, 3] / moments.higher[, 2]^1.5
  moments.kurtosis <- moments.higher[, 4] / moments.higher[, 2]^2
  moments <- list(
    mean = moments.means, var = moments.var, W = moments.W,
    B = moments.B, corr = moments.corr, Rtr = moments.Rtr,
    Rdet = moments.Rdet, skewness = moments.skewness,
    kurtosis = moments.kurtosis
  )
  return(moments)
}
