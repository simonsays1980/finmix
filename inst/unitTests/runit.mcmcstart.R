### --- Test Setup --- ###

if(TRUE) {
	## Not really needed, but can be handy 
	## when writing tests 
	library("RUnit")
	library("finmix")
}

".setUp.y" <- function() 
{
    ## Get path ##
    pkg <- "finmix"
    if (Sys.getenv("RCMDCHECK") == FALSE) {
        data.path <- file.path(getwd(), "..", 
                               "data", "poisson.data.csv")        
    } else {
        data.path <- system.file(package = pkg, 
                                 'data/poisson.data.csv')
    }
    read.csv(data.path, header = FALSE, sep = ",")
}

".setUp.S" <- function() 
{
    if (Sys.getenv("RCMDCHECK") == FALSE) {
        ind.path <- file.path(getwd(), "..", 
                              "data", 
                              "poisson.ind.csv")        
    } else {
        ind.path <- system.file(package = pkg, 
                                'data/poisson.ind.csv')
    }               
    read.csv(ind.path, header = FALSE, sep = ",")
}

## Start testing ##
"test.mcmcstart.exceptions" <- function()
{
    checkException(mcmcstart(), "check1")
    checkException(mcmcstart(fdata(), model()), "check2")
    y   <- .setUp.y()
    fdata.obj   <- fdata(y = y)
    checkException(mcmcstart(fdata.obj), "check3")
    checkException(mcmcstart(fdata.obj, model(), model()), "check4")
}

"test.mcmcstart.poisson" <- function()
{
    ## Setup
    y   <- .setUp.y()
    fdata.obj   <- fdata(y = y)
    model.obj   <- model(dist = "poisson", K = 2)
    (fdata.obj~model.obj~mcmc.obj) %=% mcmcstart(fdata.obj, model.obj)
    checkTrue(class(mcmc.obj) == "mcmc", "check1")
    checkTrue(hasS(fdata.obj), "check2")
    checkTrue(!hasPar(model.obj), "check3")
    ## Starting with indicators (@startpar = FALSE)
    fdata.obj   <- fdata(y = y)
    mcmc.obj    <- mcmc(startpar = FALSE)
    (fdata.obj~model.obj~prior.obj) %=% mcmcstart(fdata.obj, model.obj, mcmc.obj)
    checkTrue(hasPar(model.obj), "check4")
    checkTrue(!hasS(fdata.obj), "check5")
    ## Check model with fixed indicators (@indicfix = TRUE)
    setIndicfix(model.obj) <- TRUE
    fdata.obj   <- fdata(y = y)
    (fdata.obj~model.obj~mcmc.obj) %=% mcmcstart(fdata.obj, model.obj)
    checkTrue(!hasS(fdata.obj), "check6")
    options(warn = 2)
    checkException(mcmcstart(fdata.obj, model.obj), "check7")
    options(warn = 1)
}
