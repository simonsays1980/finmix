################################################################
# Copyright (c) 2013 All Rights Reserved
# author:	Barry Rowlingson
# created: January 2013
#
# avialability: https://gist.github.com/spacedman/4543212
#################################################################

#' Unstructuring assignments
#' 
#' \code{unsass} assigns multiple objects in its argument \code{rhs} (right-hand side) 
#' to multiple objects (names) chained in its argument \code{lhs} (left-hand-side).
#' 
#' This is a helper function to simplify the use of the package. The right-hand side can
#' be a function that returns multiple objects and the left-hand side must be a formula 
#' with objects (names) chained by \code{~}. Assignment works via \code{lhs %=% rhs}.
#' 
#' @param lhs A \code{formula} chaining multiple objects (names) together by \code{~}. 
#' These are the objects (names) the right-hand side should be assigned to. 
#' @param rhs A \code{list} of objects that should be assigned to the left-hand side 
#' \code{lhs}.
#' 
#' @return None. 
#' 
#' @example 
#' f_model <- model(K=2, dist= 'poisson', par=list(lambda=c(0.17, 0.12)), weight = matrix(c(0.6, 0.4), nrow=1))
#' f_data <- simulate(model)
#' mcmc <- mcmc()
#' f_data~f_model~mcmc) %=% mcmcstart(f_data, f_model, mcmc)
#' 
#' @seealso \code{mcmcstart}
"unsass" <- function(lhs, rhs) {
  nvalues <- length(rhs)
  lhss <- .getFormulaNames(lhs)
  if (length(lhss) != nvalues) {
    stop("Wrong number of values to unpack")
  }

  for (i in 1:nvalues) {
    eval(substitute(
      target <- value,
      list(target = lhss[[i]], value = rhs[[i]])
    ),
    envir = parent.frame()
    )
  }
  invisible(0)
}

assign("%=%", unsass)

".getFormulaNames" <- function(formula) {
  ## extract elements from a~b[1]~c~d
  ## recursive - might be an easier way...
  ##
  if (is.name(formula)) {
    return(formula)
  } else {
    if (is.call(formula)) {
      if (formula[[1]] == "~") {
        return(c(getFormulaNames(formula[[2]]), getFormulaNames(formula[[3]])))
      } else {
        return(formula)
      }
    }
  }
}
