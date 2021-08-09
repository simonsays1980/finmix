if(TRUE) {
	## Not really needed, but can be handy 
	## when writing tests 
	library("RUnit")
	library("finmix")
}

## Testing
"test.mcmc.default" <- function()
{
    ## Default
    mcmc.obj    <- mcmc()
    checkTrue(class(mcmc.obj) == "mcmc", "check1")
    checkTrue(mcmc.obj@startpar, "check2")
    checkTrue(mcmc.obj@ranperm, "check3")
    checkTrue(mcmc.obj@storepost, "check4")
    checkEquals(mcmc.obj@burnin, 0)
    checkEquals(mcmc.obj@M, 5000)
    checkEquals(mcmc.obj@storeS, 1000)
}

"test.mcmc.validity" <- function()
{
    ## Setup
    ## checkException 
    checkException(mcmc(burnin = -10), "check1")
    checkException(mcmc(storeS = -2), "check2")
    checkException(mcmc(M = 0), "check3")
    ## Setters
    mcmc.obj    <- mcmc()
    checkException(setBurnin(mcmc.obj) <- -10, "check4")
    checkException(setStoreS(mcmc.obj) <- -10, "check5")
    checkException(setM(mcmc.obj) <- 0, "check6")
}
