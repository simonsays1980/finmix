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

## Testing 
"test.prior.default" <- function()
{
    prior.obj <- prior()
    checkTrue(class(prior.obj) == "prior", "check1")
    checkTrue("weight" %in% slotNames(prior.obj), "check2")
    checkTrue("par" %in% slotNames(prior.obj), "check3")
    checkTrue("type" %in% slotNames(prior.obj), "check4")
    checkTrue("hier" %in% slotNames(prior.obj), "check5")
    checkTrue(all(is.na(prior.obj@weight)), "check6")
    checkEquals(length(prior.obj@par), 0)
    checkTrue(prior.obj@hier, "check7")
    checkTrue(prior.obj@type == "independent", "check8")
}

"test.prior.poisson" <- function()
{
    ## Setup
    y   <- .setUp.y()
    fdata.obj   <- fdata(y = y)
    model.obj   <- model(dist = "poisson", K = 2)
    prior.obj   <- priordefine(fdata.obj, model.obj)
    checkTrue("weight" %in% slotNames(prior.obj), "check1")
    checkTrue("type" %in% slotNames(prior.obj), "check2")
    checkTrue("par" %in% slotNames(prior.obj), "check3")
    checkTrue("hier" %in% slotNames(prior.obj), "check4")
    checkEquals(dim(prior.obj@weight), c(1, 2))
    checkEquals(length(prior.obj@par), 4)
    checkTrue(prior.obj@hier, "check5")
    checkTrue(prior.obj@type == "condconjugate", "check6")
    checkTrue(is.list(prior.obj@par), "check7")
    checkTrue("a" %in% names(prior.obj@par), "check8")
    checkTrue("b" %in% names(prior.obj@par), "check9")
    checkTrue("g" %in% names(prior.obj@par), "check10")
    checkTrue("G" %in% names(prior.obj@par), "check11")
    checkEquals(length(prior.obj@par$a), 2)
    checkEquals(length(prior.obj@par$b), 2)
    checkEquals(length(prior.obj@par$g), 1)
    checkEquals(length(prior.obj@par$G), 1)
    ## Check with no hier
    prior.var   <- prior(hier = FALSE)
    prior.obj   <- priordefine(fdata.obj, model.obj, 
                               varargin = prior.var)
    checkTrue(is.list(prior.obj@par), "check12")
    checkTrue("a" %in% names(prior.obj@par), "check13")
    checkTrue("b" %in% names(prior.obj@par), "check14")
    checkTrue(!"g" %in% names(prior.obj@par), "check15")
    checkTrue(!"G" %in% names(prior.obj@par), "check16")
    checkEquals(length(prior.obj@par$a), 2)
    checkEquals(length(prior.obj@par$b), 2)
    ## Check exceptions
    fdata.obj   <- fdata()
    checkException(priordefine(fdata.obj, model.obj), "check17")
    fdata.obj   <- fdata(y = y)
    model.obj   <- model(dist = "normult", K = 2)
    checkException(priordefine(fdata.obj, model.obj), "check18")
    model.obj   <- model(dist = "poisson", K = 2)
    checkException(priordefine(fdata.obj, model.obj, 
                               varargin = model.obj), "check19")
}
