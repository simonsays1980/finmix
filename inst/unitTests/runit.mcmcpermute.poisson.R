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

## Testing ##
"test.mcmcpermute.exceptions" <- function()
{
    mcmc.obj    <- mcmc()
    checkException(mcmcpermute(mcmc.obj), "check1")
    y   <- .setUp.y()
    S   <- .setUp.S()
    fdata.obj   <- fdata(y = y, S = S$V1)
    ## Check K == 1
    model.obj   <- model(K = 1)
    prior.obj   <- priordefine(fdata.obj, model.obj)
    mcmc.obj    <- mcmc()
    mcmcout     <- mixturemcmc(fdata.obj, model.obj,
                               prior.obj, mcmc.obj)
    options(warn = 2)
    checkException(mcmcpermute(mcmcout), "check2")
    ## Check indicfix = TRUE
    model.obj   <- model(K = 2, indicfix = TRUE)
    prior.obj   <- priordefine(fdata.obj, model.obj)
    mcmcout     <- mixturemcmc(fdata.obj, model.obj,
                               prior.obj, mcmc.obj)
    checkException(mcmcpermute(mcmcout), "check3")
    model.obj   <- model(K = 2, indicfix = FALSE)
    prior.obj   <- priordefine(fdata.obj, model.obj)
    mcmcout     <- mixturemcmc(fdata.obj, model.obj,
                               prior.obj, mcmc.obj)
    checkException(mcmcpermute(mcmcout, method = "sth"), "check3")
    checkException(mcmcpermute(mcmcout, model.obj, method = "Stephens1997b"), "check4")
    options(warn = 1)
}

"test.mcmcpermute.kmeans" <- function()
{
    ## Setup
    y   <- .setUp.y()
    S   <- .setUp.S()
    fdata.obj   <- fdata(y = y, S = S$V1)
    model.obj   <- model(K = 2)
    prior.obj   <- priordefine(fdata.obj, model.obj)
    mcmc.obj    <- mcmc()
    mcmcout     <- mixturemcmc(fdata.obj, model.obj,
                               prior.obj, mcmc.obj)
    perm        <- mcmcpermute(mcmcout)
    checkTrue(class(perm) == "mcmcoutputpermhierpost", "check1")
    checkTrue(perm@ranperm, "check2")
    checkEquals(NROW(perm@par$lambda), perm@M)
    checkEquals(NCOL(perm@par$lambda), perm@model@K)
    checkEquals(NROW(perm@log$mixlik), perm@M)
    checkEquals(NCOL(perm@log$mixlik), 1)
    checkEquals(NROW(perm@log$mixprior), perm@M)
    checkEquals(NCOL(perm@log$mixprior), 1)
    checkEquals(NROW(perm@weight), perm@M)
    checkEquals(NCOL(perm@weight), perm@model@K)
    checkEquals(NROW(perm@S), fdata.obj@N)
    checkEquals(NCOL(perm@S), mcmc.obj@storeS)
    checkEquals(NROW(perm@NK), perm@M)
    checkEquals(NCOL(perm@NK), perm@model@K)
    checkEquals(NROW(perm@clust), fdata.obj@N)
    checkEquals(NCOL(perm@clust), 1)
    checkEquals(NROW(perm@ST), perm@M)
    checkEquals(NCOL(perm@ST), 1)
    checkEquals(NROW(perm@post$par$a), perm@M)
    checkEquals(NCOL(perm@post$par$a), perm@model@K)
    checkEquals(NROW(perm@post$par$b), perm@M)
    checkEquals(NCOL(perm@post$par$b), perm@model@K)
    checkEquals(NROW(perm@post$weight), perm@M)
    checkEquals(NCOL(perm@post$weight), perm@model@K)
    checkEquals(NROW(perm@hyper$b), perm@M)    
    checkEquals(NCOL(perm@hyper$b), 1) 
    ## perm attributes
    checkTrue(perm@Mperm <= perm@M, "check3")
    checkEquals(NROW(perm@parperm$lambda), perm@Mperm)
    checkEquals(NCOL(perm@parperm$lambda), perm@model@K)
    checkEquals(NROW(perm@logperm$mixlik), perm@Mperm)
    checkEquals(NCOL(perm@logperm$mixlik), 1)
    checkEquals(NROW(perm@logperm$mixprior), perm@Mperm)
    checkEquals(NCOL(perm@logperm$mixprior), 1)
    checkEquals(NROW(perm@weightperm), perm@Mperm)
    checkEquals(NCOL(perm@weightperm), perm@model@K)
    checkEquals(NROW(perm@Sperm), fdata.obj@N)
    checkEquals(NCOL(perm@Sperm), mcmc.obj@storeS)
    checkEquals(NROW(perm@NKperm), perm@Mperm)
    checkEquals(NCOL(perm@NKperm), perm@model@K)
    checkEquals(NROW(perm@clust), fdata.obj@N)
    checkEquals(NCOL(perm@clust), 1)
    checkEquals(NROW(perm@STperm), perm@Mperm)
    checkEquals(NCOL(perm@STperm), 1)
    checkEquals(NROW(perm@postperm$par$a), perm@Mperm)
    checkEquals(NCOL(perm@postperm$par$a), perm@model@K)
    checkEquals(NROW(perm@postperm$par$b), perm@Mperm)
    checkEquals(NCOL(perm@postperm$par$b), perm@model@K)
    checkEquals(NROW(perm@postperm$weight), perm@Mperm)
    checkEquals(NCOL(perm@postperm$weight), perm@model@K)
}

"test.mcmcpermute.Stephens1997a" <- function()
{
    ## Setup
    y   <- .setUp.y()
    S   <- .setUp.S()
    fdata.obj   <- fdata(y = y, S = S$V1)
    model.obj   <- model(K = 2)
    prior.obj   <- priordefine(fdata.obj, model.obj)
    mcmc.obj    <- mcmc()
    mcmcout     <- mixturemcmc(fdata.obj, model.obj,
                               prior.obj, mcmc.obj)
    perm        <- mcmcpermute(mcmcout, method = "Stephens1997a")
    checkTrue(class(perm) == "mcmcoutputpermhierpost", "check1")
    checkTrue(perm@ranperm, "check2")
    checkEquals(NROW(perm@par$lambda), perm@M)
    checkEquals(NCOL(perm@par$lambda), perm@model@K)
    checkEquals(NROW(perm@log$mixlik), perm@M)
    checkEquals(NCOL(perm@log$mixlik), 1)
    checkEquals(NROW(perm@log$mixprior), perm@M)
    checkEquals(NCOL(perm@log$mixprior), 1)
    checkEquals(NROW(perm@weight), perm@M)
    checkEquals(NCOL(perm@weight), perm@model@K)
    checkEquals(NROW(perm@S), fdata.obj@N)
    checkEquals(NCOL(perm@S), mcmc.obj@storeS)
    checkEquals(NROW(perm@NK), perm@M)
    checkEquals(NCOL(perm@NK), perm@model@K)
    checkEquals(NROW(perm@clust), fdata.obj@N)
    checkEquals(NCOL(perm@clust), 1)
    checkEquals(NROW(perm@ST), perm@M)
    checkEquals(NCOL(perm@ST), 1)
    checkEquals(NROW(perm@post$par$a), perm@M)
    checkEquals(NCOL(perm@post$par$a), perm@model@K)
    checkEquals(NROW(perm@post$par$b), perm@M)
    checkEquals(NCOL(perm@post$par$b), perm@model@K)
    checkEquals(NROW(perm@post$weight), perm@M)
    checkEquals(NCOL(perm@post$weight), perm@model@K)
    checkEquals(NROW(perm@hyper$b), perm@M)    
    checkEquals(NCOL(perm@hyper$b), 1) 
    ## perm attributes
    checkTrue(perm@Mperm <= perm@M, "check3")
    checkEquals(NROW(perm@parperm$lambda), perm@Mperm)
    checkEquals(NCOL(perm@parperm$lambda), perm@model@K)
    checkEquals(NROW(perm@logperm$mixlik), perm@Mperm)
    checkEquals(NCOL(perm@logperm$mixlik), 1)
    checkEquals(NROW(perm@logperm$mixprior), perm@Mperm)
    checkEquals(NCOL(perm@logperm$mixprior), 1)
    checkEquals(NROW(perm@weightperm), perm@Mperm)
    checkEquals(NCOL(perm@weightperm), perm@model@K)
    checkEquals(NROW(perm@Sperm), fdata.obj@N)
    checkEquals(NCOL(perm@Sperm), mcmc.obj@storeS)
    checkEquals(NROW(perm@NKperm), perm@Mperm)
    checkEquals(NCOL(perm@NKperm), perm@model@K)
    checkEquals(NROW(perm@clust), fdata.obj@N)
    checkEquals(NCOL(perm@clust), 1)
    checkEquals(NROW(perm@STperm), perm@Mperm)
    checkEquals(NCOL(perm@STperm), 1)
    checkEquals(NROW(perm@postperm$par$a), perm@Mperm)
    checkEquals(NCOL(perm@postperm$par$a), perm@model@K)
    checkEquals(NROW(perm@postperm$par$b), perm@Mperm)
    checkEquals(NCOL(perm@postperm$par$b), perm@model@K)
    checkEquals(NROW(perm@postperm$weight), perm@Mperm)
    checkEquals(NCOL(perm@postperm$weight), perm@model@K)
}

"test.mcmcpermute.Stephens1997b" <- function()
{
    ## Setup
    y   <- .setUp.y()
    S   <- .setUp.S()
    fdata.obj   <- fdata(y = y, S = S$V1)
    model.obj   <- model(K = 2)
    prior.obj   <- priordefine(fdata.obj, model.obj)
    mcmc.obj    <- mcmc()
    mcmcout     <- mixturemcmc(fdata.obj, model.obj,
                               prior.obj, mcmc.obj)
    perm        <- mcmcpermute(mcmcout, fdata.obj, method = "Stephens1997b")
    checkTrue(class(perm) == "mcmcoutputpermhierpost", "check1")
    checkTrue(perm@ranperm, "check2")
    checkEquals(NROW(perm@par$lambda), perm@M)
    checkEquals(NCOL(perm@par$lambda), perm@model@K)
    checkEquals(NROW(perm@log$mixlik), perm@M)
    checkEquals(NCOL(perm@log$mixlik), 1)
    checkEquals(NROW(perm@log$mixprior), perm@M)
    checkEquals(NCOL(perm@log$mixprior), 1)
    checkEquals(NROW(perm@weight), perm@M)
    checkEquals(NCOL(perm@weight), perm@model@K)
    checkEquals(NROW(perm@S), fdata.obj@N)
    checkEquals(NCOL(perm@S), mcmc.obj@storeS)
    checkEquals(NROW(perm@NK), perm@M)
    checkEquals(NCOL(perm@NK), perm@model@K)
    checkEquals(NROW(perm@clust), fdata.obj@N)
    checkEquals(NCOL(perm@clust), 1)
    checkEquals(NROW(perm@ST), perm@M)
    checkEquals(NCOL(perm@ST), 1)
    checkEquals(NROW(perm@post$par$a), perm@M)
    checkEquals(NCOL(perm@post$par$a), perm@model@K)
    checkEquals(NROW(perm@post$par$b), perm@M)
    checkEquals(NCOL(perm@post$par$b), perm@model@K)
    checkEquals(NROW(perm@post$weight), perm@M)
    checkEquals(NCOL(perm@post$weight), perm@model@K)
    checkEquals(NROW(perm@hyper$b), perm@M)    
    checkEquals(NCOL(perm@hyper$b), 1) 
    ## perm attributes
    checkTrue(perm@Mperm <= perm@M, "check3")
    checkEquals(NROW(perm@parperm$lambda), perm@Mperm)
    checkEquals(NCOL(perm@parperm$lambda), perm@model@K)
    checkEquals(NROW(perm@logperm$mixlik), perm@Mperm)
    checkEquals(NCOL(perm@logperm$mixlik), 1)
    checkEquals(NROW(perm@logperm$mixprior), perm@Mperm)
    checkEquals(NCOL(perm@logperm$mixprior), 1)
    checkEquals(NROW(perm@weightperm), perm@Mperm)
    checkEquals(NCOL(perm@weightperm), perm@model@K)
    checkEquals(NROW(perm@Sperm), fdata.obj@N)
    checkEquals(NCOL(perm@Sperm), mcmc.obj@storeS)
    checkEquals(NROW(perm@NKperm), perm@Mperm)
    checkEquals(NCOL(perm@NKperm), perm@model@K)
    checkEquals(NROW(perm@clust), fdata.obj@N)
    checkEquals(NCOL(perm@clust), 1)
    checkEquals(NROW(perm@STperm), perm@Mperm)
    checkEquals(NCOL(perm@STperm), 1)
    checkEquals(NROW(perm@postperm$par$a), perm@Mperm)
    checkEquals(NCOL(perm@postperm$par$a), perm@model@K)
    checkEquals(NROW(perm@postperm$par$b), perm@Mperm)
    checkEquals(NCOL(perm@postperm$par$b), perm@model@K)
    checkEquals(NROW(perm@postperm$weight), perm@Mperm)
    checkEquals(NCOL(perm@postperm$weight), perm@model@K)
}
