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
"test.groupmoments.discrete" <- function()
{
    ## Setup
    y   <- .setUp.y() 
    S   <- .setUp.S()
    fdata.obj   <- fdata(y = y, S = S$V1)
    mom         <- groupmoments(fdata.obj)
    checkTrue(class(mom) == "groupmoments", "check1") 
    checkTrue("NK" %in% slotNames(mom), "check2")
    checkTrue("mean" %in% slotNames(mom), "check3")
    checkTrue("WK" %in% slotNames(mom), "check4")
    checkTrue("var" %in% slotNames(mom), "check5")
    checkTrue("fdata" %in% slotNames(mom), "check6")
    checkTrue(!all(is.na(mom@NK)), "check7")
    checkTrue(!all(is.na(mom@mean)), "check8")
    checkTrue(!all(is.na(mom@WK)), "check9")
    checkTrue(!all(is.na(mom@var)), "check10")
    checkEquals(dim(mom@NK), 2)
    checkEquals(dim(mom@mean), c(1, 2))
    checkEquals(dim(mom@WK), c(1, 1, 2))
    checkEquals(dim(mom@var), c(1, 1, 2))
    checkTrue(class(mom@fdata) == "fdata", "check11")
    ## Check exception
    fdata.obj   <- fdata()
    checkException(groupmoments(fdata.obj), "check12")
    fdata.obj   <- fdata(y = y) 
    checkException(groupmoments(fdata.obj), "check13")
}

"test.sdatamoments.discrete" <- function()
{
    ## Setup
    y   <- .setUp.y()
    S   <- .setUp.S()
    fdata.obj   <- fdata(y = y, S = S$V1)
    mom         <- sdatamoments(fdata.obj)
    checkTrue(class(mom) == "sdatamoments", "check1")
    checkTrue("gmoments" %in% slotNames(mom), "check2")
    checkTrue("fdata" %in% slotNames(mom), "check3")
    checkTrue(class(mom@gmoments) == "groupmoments", "check4")
    checkTrue(class(mom@fdata) == "fdata", "check5")
    checkTrue(!all(is.na(mom@gmoments@NK)), "check6")
    checkTrue(!all(is.na(mom@gmoments@mean)), "check7")
    checkTrue(!all(is.na(mom@gmoments@WK)), "check8")
    checkTrue(!all(is.na(mom@gmoments@var)), "check9")
    checkEquals(dim(mom@gmoments@NK), 2)
    checkEquals(dim(mom@gmoments@mean), c(1, 2))
    checkEquals(dim(mom@gmoments@WK), c(1, 1, 2))
    checkEquals(dim(mom@gmoments@var), c(1, 1, 2))
    ## Check exception
    fdata.obj   <- fdata()
    checkException(sdatamoments(fdata.obj), "check10")
    fdata.obj   <- fdata(y = y)
    checkException(sdatamoments(fdata.obj), "check11")
}

"test.datamoments.discrete" <- function()
{
    ## Setup
    ## No indicators
    y   <- .setUp.y()
    fdata.obj   <- fdata(y = y)
    mom <- datamoments(fdata.obj)
    checkTrue(class(mom) == "ddatamoments", "check1")
    checkTrue("mean" %in% slotNames(mom), "check2")
    checkTrue("var" %in% slotNames(mom), "check3") 
    checkTrue("fdata" %in% slotNames(mom), "check4") 
    checkTrue("factorial" %in% slotNames(mom), "check6")
    checkTrue("over" %in% slotNames(mom), "check7")
    checkTrue("zero" %in% slotNames(mom), "check8")
    checkTrue("smoments" %in% slotNames(mom), "check9") 
    checkTrue(!all(is.na(mom@mean)), "check10")
    checkTrue(!all(is.na(mom@var)), "check11")
    checkTrue(!all(is.na(mom@factorial)), "check12")
    checkTrue(!is.na(mom@over), "check13")
    checkTrue(!is.na(mom@zero), "check14")
    checkEquals(length(mom@mean), 1)
    checkEquals(dim(mom@var), c(1, 1))
    checkEquals(dim(mom@factorial), c(4, 1))
    checkEquals(length(mom@over), 1)
    checkEquals(length(mom@zero), 1)
    checkTrue(class(mom@fdata) == "fdata", "check15")
    checkTrue(class(mom@smoments) == "NULL", "check16")
    ## Check exceptions
    fdata.obj <- fdata()
    checkException(datamoments(fdata.obj), "check17")
    checkException(datamoments(), "check18")
    ## Check with y and S
    S   <- .setUp.S()
    fdata.obj   <- fdata(y = y, S = S$V1)
    mom <- datamoments(fdata.obj)
    checkTrue(class(mom@smoments) == "sdatamoments", "check19")
    checkTrue(class(mom@smoments@gmoments) == "groupmoments", "check20")
    gmom    <- mom@smoments@gmoments
    checkTrue(!all(is.na(gmom@NK)), "check21")
    checkTrue(!all(is.na(gmom@mean)), "check22")
    checkTrue(!all(is.na(gmom@WK)), "check23")
    checkTrue(!all(is.na(gmom@var)), "check24")
    checkEquals(dim(gmom@NK), 2)
    checkEquals(dim(gmom@mean), c(1, 2))
    checkEquals(dim(gmom@WK), c(1, 1, 2))
    checkEquals(dim(gmom@var), c(1, 1, 2))
    checkTrue(class(gmom@fdata) == "fdata", "check25")
}

