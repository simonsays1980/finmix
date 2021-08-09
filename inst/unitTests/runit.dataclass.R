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
"test.dataclass.exceptions" <- function()
{
    y   <- .setUp.y()
    checkException(dataclass(), "check1")
    fdata.obj   <- fdata(y = y)
    checkException(dataclass(fdata.obj), "check2")
    model.obj   <- model("poisson", K = 2)
    checkException(dataclass(model.obj, model.obj), "check3")
    checkException(dataclass(fdata.obj, model.obj), "check4")
    mcmc.obj <- mcmc(startpar = FALSE) 
    (fdata.obj~model.obj~mcmc.obj) %=% mcmcstart(fdata.obj, model.obj, mcmc.obj)
    S   <- as.matrix(sample(seq(1, 3), 100, TRUE))
    fdata.obj@S <- S
    checkException(dataclass(fdata.obj, model.obj), "check5")
    ## Model with fixed indicators 
    model.obj@indicfix <- TRUE
    checkException(dataclass(fdata.obj, model.obj), "check6")
}

"test.dataclass.poisson.default" <- function()
{
    ## Setup
    y   <- .setUp.y()
    fdata.obj   <- fdata(y = y)
    model.obj   <- model("poisson", K = 2)
    mcmc.obj    <- mcmc(startpar = FALSE)
    (fdata.obj~model.obj~mcmc.obj) %=% mcmcstart(fdata.obj, model.obj, mcmc.obj)
    datac.obj   <- dataclass(fdata.obj, model.obj)
    checkTrue(class(datac.obj) == "dataclass", "check1")
    checkTrue(!all(is.na(datac.obj@logpy)), "check2")
    checkTrue(!all(is.na(datac.obj@prob)), "check3")
    checkTrue(!is.na(datac.obj@mixlik), "check4")
    checkTrue(!is.na(datac.obj@entropy), "check5")
    checkTrue(!all(is.na(datac.obj@loglikcd)), "check6")
    ## With simulated indicators
    datac.list  <- dataclass(fdata.obj, model.obj, simS = TRUE)
    datac.obj   <- datac.list$dataclass
    S           <- datac.list$S
    checkTrue(!is.na(datac.obj@postS), "check7")
    checkTrue(!all(is.na(S)), "check8")
    ## With fixed indicators
    model.obj@indicfix  <- TRUE
    datac.obj   <- dataclass(fdata.obj, model.obj) 
    checkTrue(!all(is.na(datac.obj@logpy)), "check9")
    checkTrue(all(is.na(datac.obj@prob)), "check10")
    checkTrue(all(is.na(datac.obj@mixlik)), "check11")
    checkTrue(all(is.na(datac.obj@entropy)), "check12")
    checkTrue(all(is.na(datac.obj@postS)), "check13")
    ## With fixed indicators and simS = TRUE
    datac.obj   <- dataclass(fdata.obj, model.obj, simS = TRUE)
}
