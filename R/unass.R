#' Unstructuring assignments
#'
#' @description \code{unsass} assigns multiple objects in its argument
#' \code{rhs} (right-hand side) to multiple objects (names) chained in its
#' argument \code{lhs} (left-hand-side).
#'
#' This is a helper function to simplify the use of the package. The right-hand
#' side can be a function that returns multiple objects and the left-hand side
#' must be a formula with objects (names) chained by \code{~}. Assignment works
#' via \code{lhs %=% rhs}.
#'
#' @param lhs A \code{formula} chaining multiple objects (names) together by
#'   \code{~}. These are the objects (names) the right-hand side should be
#'   assigned to.
#' @param rhs A \code{list} of objects that should be assigned to the left-hand
#'   side \code{lhs}.
#' @name unsass
#' @rdname unsass
#' 
#' @examples 
#' f_model <- model(K=2, dist= 'poisson', par=list(lambda=c(0.17, 0.12))) 
#' f_data <- simulate(f_model) 
#' mcmc <- mcmc() 
#' (f_data~f_model~mcmc) %=% mcmcstart(f_data, f_model, mcmc)
#'
#' @author Barry Rowlingson (January, 2013)
#' @seealso 
#' * \code{\link{mcmcstart}} for generating starting parameters or indicators
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

#' Assigns `\%=\%` to the `unsass()` function
#'
#'
#' @name %=%
#' @rdname unsass
#' @aliases unsass
#' @export
#' @seealso
#' * \code{\link{unsass}}
assign("%=%", unsass)

#' Extracts elements from a formula
#' 
#' @description 
#' Used in the \link{unsass} function. 
#' 
#' @param formula A formula defining the left-hand and right-hand side of the
#'   assignment.
#' @return Elements from the formula.
#' @noRd
".getFormulaNames" <- function(formula) {
  ## extract elements from a~b[1]~c~d
  ## recursive - might be an easier way...
  ##
  if (is.name(formula)) {
    return(formula)
  } else {
    if (is.call(formula)) {
      if (formula[[1]] == "~") {
        return(c(.getFormulaNames(formula[[2]]), .getFormulaNames(formula[[3]])))
      } else {
        return(formula)
      }
    }
  }
}
