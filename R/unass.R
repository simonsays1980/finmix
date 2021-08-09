# Copyright (c) 2013 All Rights Reserved 
# author:	Barry Rowlingson
# created:	January 2013
#
# This code has been copied from 'https://gist.github.com/spacedman/4543212'
# and is used in package 'finmix' to assign several modified objects 
# to a list. 

unsass <- function(lhs,rhs)
{
  nvalues = length(rhs)
  lhss = getFormulaNames(lhs)
  if (length(lhss)!=nvalues) {
    stop("Wrong number of values to unpack")
  }
 
  for (i in 1:nvalues) {
    eval(substitute(target <- value,
                    list(target=lhss[[i]],value=rhs[[i]])),
         envir=parent.frame())
  }
  invisible(0)
}
 
assign("%=%",unsass)

getFormulaNames <- function(formula)
{
  ## extract elements from a~b[1]~c~d
  ## recursive - might be an easier way...
  ##
  if (is.name(formula)) {
    return(formula)
  } else {
    if (is.call(formula)) {
      if (formula[[1]]=="~") {
        return(c(getFormulaNames(formula[[2]]),getFormulaNames(formula[[3]])))
      } else {
        return(formula)
      }
    }
  }
}
