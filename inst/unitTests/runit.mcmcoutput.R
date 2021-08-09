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
"test.mcmcoutput.subseq" <- function()
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
    index       <- matrix(seq(1, mcmcout@M) < 6, nrow = mcmcout@M,
                          ncol = 1)
    mcmcsub     <- subseq(mcmcout, index)
    checkEquals(mcmcsub@M, 5)
    ## Test exceptions
    checkException(subseq(mcmcout), "check1")
    checkException(subseq(mcmcout, TRUE), "check2")
    index       <- matrix(1, nrow = mcmcout@M, ncol = 2)
    checkException(subseq(mcmcout, index), "check3")
    index       <- matrix(mcmcout@M > 1,  nrow = mcmcout@M, ncol = 2)
    checkException(subseq(mcmcout, index), "check4")
}

"test.mcmcoutput.swapElements" <- function()
{
    ## Setup
    y   <- .setUp.y()
    S   <- .setUp.S()
    fdata.obj   <- fdata(y = y, S = S$V1)
    model.obj   <- model(K = 2)
    prior.obj   <- priordefine(fdata.obj, model.obj)
    mcmc.obj    <- mcmc(storeS = 1)
    mcmcout     <- mixturemcmc(fdata.obj, model.obj,
                               prior.obj, mcmc.obj)
    index       <- matrix(as.integer(2), nrow = 100, ncol = 2)
    checkException(swapElements(mcmcout, index), "check1")
    index       <- matrix(as.integer(c(2, 1)), nrow = mcmcout@M, 
                          ncol = 2, byrow = TRUE)
    mcmcout.perm    <- swapElements(mcmcout, index)
    checkTrue(any(mcmcout@par$lambda != mcmcout.perm@par$lambda), "check2")
    checkTrue(any(mcmcout@post$par$a != mcmcout.perm@post$par$a), "check3")
    checkTrue(any(mcmcout@post$par$b != mcmcout.perm@post$par$b), "check4")
    checkTrue(any(mcmcout@weight != mcmcout.perm@weight), "check5")
    checkTrue(any(mcmcout@ST != mcmcout.perm@ST), "check6")
    checkTrue(any(mcmcout@S != mcmcout.perm@S), "check7")
    checkTrue(any(mcmcout@NK != mcmcout.perm@NK), "check8")
}


### --- Test [[Rcpp::export]] functions --- ###

"test.swap_cc" <- function() {
    set.seed(0)
    values          <- matrix(rnorm(20), nrow = 10, ncol = 2) 
    perm.index      <- matrix(as.integer(c(2,1)), nrow = 10, ncol = 2, byrow = TRUE)
    values.perm     <- swap_cc(values, perm.index) 
    ## Test cases ##
    checkEquals(nrow(values), nrow(values.perm))
    checkEquals(ncol(values), ncol(values.perm))
    checkTrue(!any(values == values.perm), "check3")
    ## Test exception ##
    perm.index      <- matrix(as.integer(c(2,1)), nrow = 1, ncol = 2, byrow = TRUE)
    checkException(swap_cc(values, perm.index), silent = TRUE)
    ## --- Check for K = 3 --- ##
    values          <- matrix(rnorm(30), nrow = 10, ncol = 3)
    perm.index      <- matrix(as.integer(c(2, 3, 1)), nrow = 10, ncol = 3, byrow = TRUE)
    values.perm     <- swap_cc(values, perm.index)
    ## Test cases ##
    checkEquals(nrow(values), nrow(values.perm))
    checkEquals(ncol(values), ncol(values.perm))
    checkTrue(!any(values == values.perm), "check7")
}

"test.swapInteger_cc" <- function() {
    set.seed(0)
    values          <- matrix(as.integer(rpois(20, 2)), nrow = 10, ncol = 2) 
    perm.index      <- matrix(as.integer(c(2,1)), nrow = 10, ncol = 2, byrow = TRUE)
    values.perm     <- swapInteger_cc(values, perm.index) 
    ## Test cases ##
    checkEquals(nrow(values), nrow(values.perm))
    checkEquals(ncol(values), ncol(values.perm))
    checkTrue(!any(values == values.perm), "check3")
    ## Test exception ##
    perm.index      <- matrix(as.integer(c(2,1)), nrow = 1, ncol = 2, byrow = TRUE)
    checkException(swapInteger_cc(values, perm.index), silent = TRUE)
    ## --- Check for K = 3 --- ##
    values          <- matrix(as.integer(rpois(30, 2)), nrow = 10, ncol = 3)
    perm.index      <- matrix(as.integer(c(2, 3, 1)), nrow = 10, ncol = 3, byrow = TRUE)
    values.perm     <- swapInteger_cc(values, perm.index)
    ## Test cases ##
    checkEquals(nrow(values), nrow(values.perm))
    checkEquals(ncol(values), ncol(values.perm))
    checkTrue(!all(values == values.perm), "check7")
}

"test.swapInd_cc" <- function() {
    set.seed(0) 
    indicator       <- matrix(sample(c(1,2), 10, replace = TRUE))
    perm.index      <- matrix(as.integer(c(2,1)), nrow = 1, ncol = 2, byrow = TRUE)
    indicator.perm  <- swapInd_cc(indicator, perm.index)
    ## Test cases ##
    checkEquals(nrow(indicator.perm), nrow(indicator))
    checkEquals(ncol(indicator.perm), ncol(indicator))
    checkTrue(!any(indicator == indicator.perm), "check3")
    ## Test exception ##
    perm.index      <- matrix(c(2,1), nrow = 2, ncol = 2, byrow = TRUE)
    checkException(swapInd_cc(indicator, perm.index), silent = TRUE)
    ## --- Check K = 3 --- ##
    set.seed(0)
    indicator       <- matrix(sample(c(1, 2, 3), 10, replace = TRUE))
    perm.index      <- matrix(as.integer(c(2, 3, 1)), nrow = 1, ncol = 3, byrow = TRUE)
    indicator.perm  <- swapInd_cc(indicator, perm.index) 
    ## Test cases ##
    checkEquals(nrow(indicator), nrow(indicator.perm))
    checkEquals(ncol(indicator), ncol(indicator.perm))
    checkTrue(!any(indicator == indicator.perm), "check7")
}

"test.swapST_cc" <- function() {
    set.seed(0)
    indicator       <- matrix(sample(c(1,2), 10, replace = TRUE))
    perm.index      <- matrix(as.integer(c(2,1)), nrow = 10, ncol = 2, byrow = TRUE)
    indicator.perm  <- swapST_cc(indicator, perm.index)
    ## Test cases ##
    checkEquals(nrow(indicator), nrow(indicator.perm))
    checkEquals(ncol(indicator), ncol(indicator.perm))
    checkTrue(!any(indicator == indicator.perm), "check3")
    ## Test exception ##
    perm.index      <- matrix(as.integer(c(2,1)), nrow = 2, ncol = 2, byrow = TRUE)
    checkException(swapST_cc(indicator, perm.index), silent = TRUE)
    ## --- Check for K = 3 --- ##
    set.seed(0) 
    indicator       <- matrix(sample(c(1, 2, 3), 10, replace = TRUE))
    perm.index      <- matrix(as.integer(c(2, 3, 1)), nrow = 10, ncol = 3, byrow = TRUE)
    indicator.perm  <- swapST_cc(indicator, perm.index)
    ## Test cases ##
    checkEquals(nrow(indicator), nrow(indicator.perm))
    checkEquals(ncol(indicator), ncol(indicator.perm))
    checkTrue(!any(indicator == indicator.perm), "check7")
}

