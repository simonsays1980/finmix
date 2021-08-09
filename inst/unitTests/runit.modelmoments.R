### --- Test Setup --- ###

if(TRUE) {
	## Not really needed, but can be handy 
	## when writing tests 
	library("RUnit")
	library("finmix")
}

## Start Testing ##

"test.modelmoments.poisson" <- function()
{
    ## Default 
    model.obj   <- model()
    checkException(mom <- modelmoments(model.obj), "check1")
    pars    <- list(lambda = 112)
    setPar(model.obj) <- pars
    mom <- modelmoments(model.obj)
    checkTrue("mean" %in% slotNames(mom), "check2")
    checkTrue("var" %in% slotNames(mom), "check3")
    checkTrue("factorial" %in% slotNames(mom), "check4")
    checkTrue("over" %in% slotNames(mom), "check5")
    checkTrue("zero" %in% slotNames(mom), "check6")
    checkEquals(NROW(mom@mean), 1)
    checkEquals(NCOL(mom@mean), 1)
    checkEquals(nrow(mom@var), 1)
    checkEquals(ncol(mom@var), 1)
    checkEquals(dim(mom@factorial)[1], 4)
    checkEquals(dim(mom@factorial)[2], 1)
    checkEquals(NROW(mom@over), 1)
    checkEquals(NCOL(mom@over), 1)
    checkEquals(NROW(mom@zero), 1)
    checkEquals(NCOL(mom@zero), 1)
    ## Check K = 2
    setK(model.obj) <- 2
    setPar(model.obj) <- list(lambda = c(4, 5))
    mom     <- modelmoments(model.obj)
    setWeight(model.obj)    <- matrix()
    checkException(modelmoments(model.obj), "check2")
}

