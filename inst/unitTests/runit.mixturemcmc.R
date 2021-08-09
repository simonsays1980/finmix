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

## Testing
"test.mixturemcmc.exception" <- function()
{
    ## Default
    checkException(mixturemcmc(), "check1")
    y   <- .setUp.y()
    fdata.obj   <- fdata(y = y)
    checkException(mixturemcmc(fdata.obj, fdata.obj, 
                               fdata.obj, fdata.obj), "check2")
    model.obj   <- model(dist = "poisson", K = 2)
    checkException(mixturemcmc(fdata.obj, model.obj, 
                               model.obj, model.obj), "check3") 
    prior.obj   <- priordefine(fdata.obj, model.obj)
    checkException(mixturemcmc(fdata.obj, model.obj,
                               prior.obj, prior.obj), "check4")
    ## Check empty data slot.
    fdata.obj   <- fdata() 
    mcmc.obj    <- mcmc()
    checkException(mixturemcmc(fdata.obj, model.obj, 
                               prior.obj, mcmc.obj), "check5")
    ## Check startpar = TRUE without starting indicators.
    fdata.obj   <- fdata(y = y)
    checkException(mixturemcmc(fdata.obj, model.obj, 
                               prior.obj, mcmc.obj), "check6")
    ## Check startpar = FALSE without starting parameters.
    mcmc.obj@startpar   <- FALSE
    checkException(mixturemcmc(fdata.obj, model.obj, 
                               prior.obj, mcmc.obj), "check7")
    ## Check fixed indicator model without indicators
    model.obj@indicfix  <- TRUE
    checkException(mixturemcmc(fdata.obj, model.obj, 
                               prior.obj, mcmc.obj), "check8")
    ## Check fixed indicator model without indicators but startpar set to TRUE
    mcmc.obj@startpar   <- TRUE
    checkException(mixturemcmc(fdata.obj, model.obj,
                               prior.obj, mcmc.obj), "check9")
    ## Check without prior parameters 
    model.obj@indicfix  <- FALSE
    prior.obj   <- prior()
    S           <- .setUp.S()
    fdata.obj   <- fdata(y = y, S = S$V1)
    checkException(mixturemcmc(fdata.obj, model.obj, 
                               prior.obj, mcmc.obj), "check10")
    ## Check wihtout prior weight 
    prior.obj   <- priordefine(fdata.obj, model.obj)
    prior.obj@weight    <- matrix()
    checkException(mixturemcmc(fdata.obj, model.obj,
                               prior.obj, mcmc.obj), "check11")
}
