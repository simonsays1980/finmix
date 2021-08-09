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
"test.fdata" <- function() 
{
    ## Default ##
    fdata.obj   <- fdata()
    checkTrue(all(is.na(fdata.obj@y)), "check1")
    checkTrue(all(is.na(fdata.obj@S)), "check2")
    checkTrue(all(is.na(fdata.obj@exp)),"check3")
    checkTrue(all(is.na(fdata.obj@T)), "check4")
    checkTrue(fdata.obj@bycolumn, "check5")
    checkTrue(!fdata.obj@sim, "check6")
    checkEquals(fdata.obj@N, 1)
    checkEquals(fdata.obj@r, 1)
    checkEquals(fdata.obj@type, "discrete")
    checkEquals(length(fdata.obj@name), 0)
}

"test.fdata.check.y" <- function()
{
    ## Setup ##
    y           <- .setUp.y()
    fdata.obj   <- fdata(y = y)
    checkTrue(!all(is.na(fdata.obj@y)), "check1")
    checkEquals(fdata.obj@N, nrow(fdata.obj@y))
    checkEquals(fdata.obj@r, ncol(fdata.obj@y))
    checkTrue(fdata.obj@bycolumn, "check2")
    ## Check row-ordering ##
    y           <- t(.setUp.y())
    fdata.obj   <- fdata(y = y)
    checkTrue(!all(is.na(fdata.obj@y)), "check3")
    checkEquals(fdata.obj@N, ncol(fdata.obj@y))
    checkEquals(fdata.obj@r, nrow(fdata.obj@y))
    checkTrue(!fdata.obj@bycolumn, "check4")
    ## Check exception
    y   <- matrix("", nrow = 20)
    checkException(fdata(y = y), "check5")
}

"test.fdata.check.N" <- function()
{
    ## Setup
    fdata.obj   <- fdata(N = 200)
    checkEquals(fdata.obj@N, 200)
    y           <- .setUp.y()
    fdata.obj   <- fdata(y = y, N = 100) 
    ## Check exception 
    checkException(fdata(y = y, N = 200), "check1")
    ## Check row-ordering
    y <- t(y)
    fdata.obj   <- fdata(y = y, N = 100)
    checkEquals(fdata.obj@N, 100)
    checkException(fdata(y = y, N = 200), "check2")
}

"test.fdata.check.r" <- function()
{
    ## Setup
    fdata.obj   <- fdata(r = 2, type = "continuous")
    checkEquals(fdata.obj@r, 2)
    y           <- .setUp.y()
    fdata.obj   <- fdata(y = y, r = 1) 
    ## Check exception 
    checkException(fdata(y = y, r = 2), "check1")
    ## Check row-ordering
    y <- t(y)
    fdata.obj   <- fdata(y = y, r = 1)
    checkEquals(fdata.obj@r, 1)
    checkException(fdata(y = y, r = 2), "check2")
}

"test.fdata.check.type" <- function()
{
    checkException(fdata(type = "jump"), "check1")
}

"test.fdata.check.S" <- function()
{
    S   <- .setUp.S()
    fdata.obj   <- fdata(S = S$V1)
    checkTrue(!all(is.na(fdata.obj@S)), "check1")
    checkEquals(fdata.obj@N, NROW(S$V1))
    checkEquals(fdata.obj@r, 1)
    ## Check row-ordering
    S <- t(S$V1)
    fdata.obj <- fdata(S = S)
    checkTrue(!all(is.na(fdata.obj@S)), "check2")
    checkEquals(fdata.obj@N, NCOL(S))
    checkEquals(fdata.obj@r, 1)
    ## Check Exception
    S   <- matrix("", nrow = 10)
    checkException(fdata(S = S), "check23")
    S   <- matrix(c(2.3, 4.1, 2.3))
    fdata.obj <- fdata(S = S)
    checkEquals(fdata.obj@S[1], 2)
    S   <- .setUp.S()
    S   <- S$V1[1:50]
    y   <- .setUp.y()
    checkException(fdata(y = y, S = S), "check4")
    S   <- .setUp.S()
    S   <- cbind(S$V1, S$V1)
    checkException(fdata(S = S), "check5")
    checkException(fdata(S = t(S)), "check6")
    S   <- c(2, 1, 1, 1, 2, -1)
    checkException(fdata(S = S), "check7")
}

"test.fdata.check.T" <- function()
{
    T   <- .setUp.S()
    fdata.obj   <- fdata(T = T$V1)
    checkTrue(!all(is.na(fdata.obj@T)), "check1")
    checkEquals(fdata.obj@N, NROW(T))
    checkEquals(fdata.obj@r, 1)
    ## Check row-ordering 
    T   <- t(T$V1)
    fdata.obj   <- fdata(T = T)
    checkTrue(!all(is.na(fdata.obj@T)), "check2")
    checkEquals(fdata.obj@N, NCOL(T))
    checkEquals(fdata.obj@r, 1)
    ## Check exceptions
    T   <- matrix("", nrow = 10)
    checkException(fdata(T = T), "check3")
    T   <- matrix(c(2.3, 4.1, 2.3))
    fdata.obj   <- fdata(T = T)    
    checkEquals(fdata.obj@T[1], 2)
    T   <- .setUp.S()
    T   <- T$V1[1:50]
    y   <- .setUp.y()
    checkException(fdata(y = y, T = T), "check4")
    T   <- .setUp.S()
    T   <- cbind(T$V1, T$V1)
    checkException(fdata(T = T), "check5")
    checkException(fdata(T = t(T)), "check6")
    T   <- c(2, 1, 2, 2, 0)
    checkException(fdata(T = T), "check7")
}

"test.fdata.check.exp" <- function()
{
    expos   <- .setUp.y()
    fdata.obj   <- fdata(exp = expos$V1)
    checkTrue(!all(is.na(fdata.obj@exp)), "check1")
    checkEquals(fdata.obj@N, NROW(expos))
    checkEquals(fdata.obj@r, 1)
    ## Check row-ordering 
    expos   <- t(expos$V1)
    fdata.obj   <- fdata(exp = expos)
    checkTrue(!all(is.na(fdata.obj@exp)), "check2")
    checkEquals(fdata.obj@N, NCOL(expos))
    checkEquals(fdata.obj@r, 1)
    ## Check exceptions
    expos   <- matrix("", nrow = 10)
    checkException(fdata(exp = expos), "check3")
    expos   <- .setUp.y()
    expos   <- expos$V1[1:50]
    y   <- .setUp.y()
    checkException(fdata(y = y, exp = expos), "check4")
    expos   <- .setUp.y()
    expos   <- cbind(expos$V1, expos$V1) 
    checkException(fdata(exp = expos), "check5")
    checkException(fdata(exp = t(expos)), "check6")
    expos   <- c(2, -1, 3, 1, 2, 0.0003)
    checkException(fdata(exp = expos), "check7")
}

"test.fdata.setY" <- function()
{
    ## Default
    fdata.obj <- fdata()
    y   <- .setUp.y()
    setY(fdata.obj) <- y
    checkTrue(!all(is.na(fdata.obj@y)), "check1")
    checkEquals(fdata.obj@N, NROW(y))
    checkEquals(fdata.obj@r, 1)
    ## Check row-ordering
    setY(fdata.obj) <- t(y)
    ## Check with S
    S   <- .setUp.S()   
    fdata.obj <- fdata(S = S$V1)
    checkEquals(fdata.obj@N, NROW(S$V1))
    setY(fdata.obj) <- y
    checkTrue(!all(is.na(fdata.obj@y)), "check2")
    checkEquals(fdata.obj@N, NROW(S$V1))
    checkEquals(fdata.obj@r, 1)
    setY(fdata.obj) <- t(y)
    checkEquals(nrow(fdata.obj@y), NROW(S$V1))
    checkEquals(ncol(fdata.obj@y), 1)
    y   <- cbind(y, y)
    setType(fdata.obj) <- "continuous"
    setY(fdata.obj) <- y
    setY(fdata.obj) <- t(y)
    ## Check exception
    y   <- matrix("", nrow = 10)
    checkException(setY(fdata.obj) <- y, "check3")
}

"test.fdata.setBycolumn" <- function()
{
    ## Default
    fdata.obj <- fdata()
    setBycolumn(fdata.obj) <- FALSE
    y   <- .setUp.y()
    fdata.obj <- fdata(y = y)
    setBycolumn(fdata.obj) <- TRUE
    checkTrue(getBycolumn(fdata.obj), "check1")
    setBycolumn(fdata.obj) <- FALSE
    checkTrue(!getBycolumn(fdata.obj), "check2")
    checkEquals(nrow(fdata.obj@y), NCOL(y))
    checkEquals(ncol(fdata.obj@y), NROW(y))
    checkEquals(fdata.obj@N, NROW(y))
    checkEquals(fdata.obj@r, NCOL(y))
    S   <- .setUp.S()
    fdata.obj <- fdata(S = S$V1)
    setBycolumn(fdata.obj) <- FALSE
    checkEquals(nrow(fdata.obj@S), NCOL(S$V1))
    checkEquals(ncol(fdata.obj@S), NROW(S$V1))
}

"test.fdata.setS" <- function()
{
    ## Default
    fdata.obj <- fdata()
    S   <- .setUp.S()
    setS(fdata.obj) <- S$V1
    checkTrue(!all(is.na(fdata.obj@S)), "check1")
    ## Check row-ordering
    setS(fdata.obj) <- t(S$V1)
    checkEquals(nrow(fdata.obj@S), NROW(S$V1))
    checkEquals(ncol(fdata.obj@S), NCOL(S$V1))
    ## Check with y
    y   <- .setUp.y()
    fdata.obj <- fdata(y = y)
    setS(fdata.obj) <- S$V1
    checkEquals(nrow(fdata.obj@S), NROW(S$V1))
    checkEquals(ncol(fdata.obj@S), NCOL(S$V1))
    ## Check with y and row-ordering
    fdata.obj <- fdata(y = t(y))
    setS(fdata.obj) <- S$V1
    checkEquals(nrow(fdata.obj@S), NCOL(S$V1))
    checkEquals(ncol(fdata.obj@S), NROW(S$V1))
    fdata.obj <- fdata(y = y)
    setS(fdata.obj) <- t(S$V1)
    checkEquals(nrow(fdata.obj@S), NROW(S$V1))
    checkEquals(ncol(fdata.obj@S), NCOL(S$V1))
    ## Check exception
    S   <- c(2, 1, 2, - 1)
    checkException(setS(fdata.obj) <- S, "check2")
    S   <- matrix("", nrow = 10)
    checkException(setS(fdata.obj) <- S, "check3")
}

"test.fdata.setExp" <- function()
{
    ## Default
    fdata.obj <- fdata()
    expos   <- .setUp.y()
    expos <- matrix(expos$V1)
    setExp(fdata.obj) <- expos
    checkTrue(!all(is.na(fdata.obj@exp)), "check1")
    ## Check row-ordering
    setExp(fdata.obj) <- t(expos)
    checkEquals(nrow(fdata.obj@exp), NROW(expos))
    checkEquals(ncol(fdata.obj@exp), NCOL(expos))
    ## Check with y
    y   <- .setUp.y()
    fdata.obj <- fdata(y = y)
    setExp(fdata.obj) <- expos
    checkEquals(nrow(fdata.obj@exp), NROW(expos))
    checkEquals(ncol(fdata.obj@exp), NCOL(expos))
    ## Check with y and row-ordering
    fdata.obj <- fdata(y = t(y))
    setExp(fdata.obj) <- expos
    checkEquals(nrow(fdata.obj@exp), NCOL(expos))
    checkEquals(ncol(fdata.obj@exp), NROW(expos))
    fdata.obj <- fdata(y = y)
    setExp(fdata.obj) <- t(expos)
    checkEquals(nrow(fdata.obj@exp), NROW(expos))
    checkEquals(ncol(fdata.obj@exp), NCOL(expos))
    ## Check exception
    expos   <- c(2, 1, 2, - 1)
    checkException(setExp(fdata.obj) <- expos, "check2")
    expos   <- matrix("", nrow = 10)
    checkException(setExp(fdata.obj) <- expos, "check3")
}

"test.fdata.setT" <- function()
{
    ## Default
    fdata.obj <- fdata()
    T   <- .setUp.S()
    setT(fdata.obj) <- T$V1
    checkTrue(!all(is.na(fdata.obj@T)), "check1")
    ## Check row-ordering
    setT(fdata.obj) <- t(T$V1)
    checkEquals(nrow(fdata.obj@T), NROW(T$V1))
    checkEquals(ncol(fdata.obj@T), NCOL(T$V1))
    ## Check with y
    y   <- .setUp.y()
    fdata.obj <- fdata(y = y)
    setT(fdata.obj) <- T$V1
    checkEquals(nrow(fdata.obj@T), NROW(T$V1))
    checkEquals(ncol(fdata.obj@T), NCOL(T$V1))
    ## Check with y and row-ordering
    fdata.obj <- fdata(y = t(y))
    setT(fdata.obj) <- T$V1
    checkEquals(nrow(fdata.obj@T), NCOL(T$V1))
    checkEquals(ncol(fdata.obj@T), NROW(T$V1))
    fdata.obj <- fdata(y = y)
    setT(fdata.obj) <- t(T$V1)
    checkEquals(nrow(fdata.obj@T), NROW(T$V1))
    checkEquals(ncol(fdata.obj@T), NCOL(T$V1))
    ## Check exception
    T   <- c(2, 1, 2, - 1)
    checkException(setT(fdata.obj) <- T, "check2")
    T   <- matrix("", nrow = 10)
    checkException(setT(fdata.obj) <- T, "check3")
}

"test.fdata.hasY" <- function()
{
    ## Default
    fdata.obj <- fdata()
    checkTrue(!hasY(fdata.obj), "check1")
    checkException(hasY(fdata.obj, verbose = TRUE), "check2")
    y   <- .setUp.y()
    fdata.obj <- fdata(y = y)
    checkTrue(hasY(fdata.obj), "check3")
}

"test.fdata.hasS" <- function()
{
    ## Default
    fdata.obj <- fdata()
    checkTrue(!hasS(fdata.obj), "check1")
    checkException(hasS(fdata.obj, verbose = TRUE), "check2")
    S   <- .setUp.S()
    fdata.obj <- fdata(S = S$V1)
    checkTrue(hasS(fdata.obj), "check3")
}

"test.fdata.hasExp" <- function()
{
    ## Default
    fdata.obj <- fdata()
    checkTrue(!hasExp(fdata.obj), "check1")
    checkException(hasExp(fdata.obj, verbose = TRUE), "check2")
    expos   <- .setUp.y()
    fdata.obj <- fdata(exp = expos$V1)
    checkTrue(hasExp(fdata.obj), "check3")
}

"test.fdata.hasT" <- function()
{
    ## Default
    fdata.obj <- fdata()
    checkTrue(!hasT(fdata.obj), "check1")
    checkException(hasT(fdata.obj, verbose = TRUE), "check2")
    T   <- .setUp.S()
    fdata.obj <- fdata(T = T$V1)
    checkTrue(hasT(fdata.obj), "check3")
}

"test.fdata.getColY" <- function()
{
    ## Default
    y   <- .setUp.y()
    fdata.obj <- fdata(y = y)
    y.out <- getColY(fdata.obj)
    checkEquals(nrow(y.out), NROW(y))
    checkEquals(ncol(y.out), NCOL(y))
}

"test.fdata.getRowY" <- function()
{
    ## Default
    y   <- .setUp.y()
    fdata.obj <- fdata(y = y)
    y.out <- getRowY(fdata.obj)
    checkEquals(ncol(y.out), NROW(y))
    checkEquals(nrow(y.out), NCOL(y))
}

"test.fdata.getColS" <- function()
{
    ## Default
    S   <- .setUp.S()
    fdata.obj <- fdata(S = S$V1)
    S.out <- getColS(fdata.obj)
    checkEquals(nrow(S.out), NROW(S))
    checkEquals(ncol(S.out), NCOL(S))
}

"test.fdata.getRowS" <- function()
{
    ## Default
    S   <- .setUp.S()
    fdata.obj <- fdata(S = S$V1)
    S.out <- getRowS(fdata.obj)
    checkEquals(ncol(S.out), NROW(S))
    checkEquals(nrow(S.out), NCOL(S))
}

"test.fdata.getColExp" <- function()
{
    ## Default
    expos   <- .setUp.y()
    fdata.obj <- fdata(exp = expos$V1)
    exp.out <- getColExp(fdata.obj)
    checkEquals(nrow(exp.out), NROW(expos))
    checkEquals(ncol(exp.out), NCOL(expos))
}

"test.fdata.getRowY" <- function()
{
    ## Default
    expos   <- .setUp.y()
    fdata.obj <- fdata(exp = expos$V1)
    exp.out <- getRowExp(fdata.obj)
    checkEquals(ncol(exp.out), NROW(expos))
    checkEquals(nrow(exp.out), NCOL(expos))
}

"test.fdata.getColT" <- function()
{
    ## Default
    T   <- .setUp.S()
    fdata.obj <- fdata(T = T$V1)
    T.out <- getColT(fdata.obj)
    checkEquals(nrow(T.out), NROW(T))
    checkEquals(ncol(T.out), NCOL(T))
}

"test.fdata.getRowT" <- function()
{
    ## Default
    T   <- .setUp.S()
    fdata.obj <- fdata(T = T$V1)
    T.out <- getRowT(fdata.obj)
    checkEquals(ncol(T.out), NROW(T))
    checkEquals(nrow(T.out), NCOL(T))
}
