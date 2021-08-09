### --- Test Setup --- ###

if(TRUE) {
	## Not really needed, but can be handy 
	## when writing tests 
	library("RUnit")
	library("finmix")
}

## Testing ##
"test.model" <- function()
{
    ## Check default ##
    model.obj <- model()
    checkTrue(all(is.na(model.obj@T)), "check2")
    checkTrue(!model.obj@indicfix, "check3")
    checkEquals(model.obj@K, 1)
    checkEquals(model.obj@r, 1)
    checkEquals(model.obj@dist, "poisson")
    checkEquals(model.obj@indicmod, "multinomial")
}

"test.model.check.dist" <- function()
{
    ## Check Exception
    checkException(model(dist = "Gamma"), "check1")
    checkException(model(indicmod = "binomial"), "check2")
}

"test.model.check.K" <- function() 
{
    ## Check Exception ## 
    checkException(model(K = -1), "check1")
    ## Check model with K = 1
    model.obj <- model(K = 1)
    checkTrue(all(is.na(model.obj@weight)), "check2")
    ## Check model with K = 2
    model.obj <- model(K = 2)
    checkTrue(!all(is.na(model.obj@weight)), "check3")
    ## Check model with K = 2 and weight
    wgt <- matrix(1/3, nrow = 1, ncol = 3)
    checkException(model(K = 2, weight = wgt), "check4")
}

"test.model.check.r" <- function()
{
    ## Check exception 
    checkException(model(r = 2))
    checkException(model(dist = "binomial", r = 2), "check1")
    checkException(model(dist = "normult", r = 1), "check2")
}

"test.model.valid.T" <- function()
{
    ## Check exception
    reps <- matrix("", nrow = 10, ncol = 1)
    checkException(model(T = reps), "check1")
    reps <- matrix(0, nrow = 10, ncol = 1)
    checkException(model(T = reps), "check2")
}

"test.model.valid.weight" <- function() 
{
    ## Check exception 
    wgt <- matrix("", nrow = 1, ncol = 2)
    checkException(model(K = 2, weight = wgt), "check1")
    ## Check K = 1 
    wgt <- matrix(1/2, nrow = 1, ncol = 2)
    checkException(model(K = 1, weight = wgt), "check2")
    storage.mode(wgt) <- "integer"
    checkException(model(weight = wgt), "check3")
}

"test.model.valid.par.poisson" <- function()
{
    ## Check exceptions
    pars <- list(mu = c(1, 2))
    checkException(model(par = pars), "check1")
    ## Check K = 2
    pars <- list(lambda = 1)
    checkException(model(K = 2, par = pars), "check2")
    ## Check with weight 
    pars <- list(lambda = c(4, 5))
    wgt <- matrix(1/3, nrow = 1, ncol = 3)
    checkException(model(par = pars, weight = wgt), "check3")
    ## Check non-numeric
    pars <- list(lambda = "")
    checkException(model(par = pars), "check4")
    ## Check negative 
    pars <- list(lambda = c(-2, 1))
    checkException(model(K = 2, par = pars), "check5")
}

## Setters 
"test.model.setDist" <- function()
{
    ## Default
    model.obj <- model() 
    checkException(setDist(model.obj) <- "Gamma", "check1")
}

"test.model.setK" <- function() 
{
    ## Default 
    model.obj <- model()
    checkException(setK(model.obj) <- -1, "check1")
}

"test.model.setWeight" <- function()
{
    ## Default 
    model.obj <- model()
    wgt <- matrix(1/2, nrow = 1, ncol = 2)
    setWeight(model.obj) <- wgt
    checkTrue(hasWeight(model.obj), "check1")
    ## NA    
    setWeight(model.obj) <- matrix()
    checkTrue(all(is.na(model.obj@weight)), "check2")
    ## Exception
    wgt <- matrix("")
    checkException(setWeight(model.obj) <- wgt, "check3")
}

"test.model.setPar" <- function()
{
    ## Default ##
    model.obj <- model()
    pars <- list(mu = 1)
    setPar(model.obj) <- pars
    checkEquals(length(model.obj@par), 1)
    checkTrue("mu" %in% names(model.obj@par), "check1")
    checkEquals(model.obj@par$mu, 1)
    ## Excpetions (turn warnings into errors)
    model.obj <- model()
    options(warn = 2)
    checkException(setPar(model.obj) <- pars, "check2")
    pars <- list(lambda = -1) 
    checkException(setPar(model.obj) <- pars, "check3")
    setK(model.obj) <- 2
    pars <- list(lambda = c(2, 3))
    setPar(model.obj) <- pars
    options(warn = 1)
}

"test.model.setT" <- function() 
{
    ## Default
    model.obj <- model()
    reps <- as.integer(c(4, 5, 1))
    setT(model.obj) <- reps
    checkTrue(!all(is.na(model.obj@T)), "check1")
    checkEquals(model.obj@T[1], 4)
    ## Exception
    reps <- matrix("")
    checkException(setT(model.obj) <- reps, "check2")
}

"test.model.hasweight" <- function() 
{
    ## Default ##
    model.obj <- model()
    checkTrue(!hasWeight(model.obj), "check1")
    checkException(hasWeight(model.obj, verbose = TRUE), "check2")
    ## Check K = 2
    model.obj <- model(K = 2)
    checkTrue(hasWeight(model.obj), "check3")
}

"test.model.haspar" <- function() 
{
    ## Default
    model.obj <- model()
    checkTrue(!hasPar(model.obj), "check1")
    setK(model.obj) <- 2
    pars <- list(lambda = c(4, 5))
    setPar(model.obj) <- pars
    checkTrue(hasPar(model.obj), "check2")
    ## Exception
    model.obj <- model()
    checkException(hasPar(model.obj, verbose = TRUE), "check3")
    pars <- list(mu = 1)
    setPar(model.obj) <- pars
    checkException(hasPar(model.obj, verbose = TRUE), "check4")
}

"test.model.hasT" <- function()
{
    ## default
    model.obj <- model()
    checkTrue(!hasT(model.obj), "check1")
    reps <- matrix(2, nrow = 1, ncol = 10)
    storage.mode(reps) <- "integer"
    setT(model.obj) <- reps
    checkTrue(hasT(model.obj), "check2")
}

## implemented later when normult model is implemented 
## test.model.mixturemar ##
