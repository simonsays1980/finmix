## Copyright (C) 2013 Lars Simon Zehnder
#
# This file is part of finmix.
#
# finmix is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# finmix is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Rcpp.  If not, see <http://www.gnu.org/licenses/>.

## Load the dynamic library
#' @useDynLib finmix
#' @importFrom Rcpp sourceCpp
NULL

## Class 'model' --------------------------------------------------
#' Simulates data from a finite mixture model
#' @export
#' @docType methods
#' @keywords internal
setGeneric("simulate", function(model, N = 100, varargin, seed = 0) standardGeneric("simulate"))

#' Plots the point process of a finite mixture model
#' @export
#' @docType methods
#' @keywords internal
#' @rdname plotPointProc-generic
setGeneric("plotPointProc", function(x, dev = TRUE, ...) standardGeneric("plotPointProc"))

#' Checks a finite mixture model for the weight
#' @export
#' @docType methods
#' @keywords internal
setGeneric("hasWeight", function(object, verbose = FALSE) standardGeneric("hasWeight"))

#' Checks a finite mixture model for repetitions
#' @export
#' @docType methods
#' @keywords internal
setGeneric("hasT", function(object, verbose = FALSE) standardGeneric("hasT"))

#' Checks a finite mixture model for the parameters
#' @export
#' @docType methods
#' @keywords internal
setGeneric("hasPar", function(object, verbose = FALSE) standardGeneric("hasPar"))

#' Extracts the marginal distribution from a finite mixture model
#' @export
#' @docType methods
#' @keywords internal 
setGeneric("mixturemar", function(object, J) standardGeneric("mixturemar"))

#' Getter for the `dist` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getDist", function(object) standardGeneric("getDist"))

#' Getter for the `r` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getR", function(object) standardGeneric("getR"))

#' Getter for the `K` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getK", function(object) standardGeneric("getK"))

#' Getter for the `weight` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getWeight", function(object) standardGeneric("getWeight"))

#' Getter for the `par` slot
#' @export 
#' @docType methods
#' @keywords internal
setGeneric("getPar", function(object) standardGeneric("getPar"))

#' Getter for the `indicmod` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getIndicmod", function(object) standardGeneric("getIndicmod"))

#' Getter for the `indicfix` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getIndicfix", function(object) standardGeneric("getIndicfix"))

#' Getter for the `T` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getT", function(object) standardGeneric("getT"))

#' Setter for the `dist` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("setDist<-", function(object, value) standardGeneric("setDist<-"))

#' Setter for the `r` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("setR<-", function(object, value) standardGeneric("setR<-"))

#' Setter for the `K` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("setK<-", function(object, value) standardGeneric("setK<-"))

#' Setter for the `weight` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("setWeight<-", function(object, value) standardGeneric("setWeight<-"))

#' Setter for the `par` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("setPar<-", function(object, value) standardGeneric("setPar<-"))

#' Setter for the `indicmod` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("setIndicmod<-", function(object, value) standardGeneric("setIndicmod<-"))

#' Setter for the `indicfix` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("setIndicfix<-", function(object, value) standardGeneric("setIndicfix<-"))

#' Setter for the `T` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("setT<-", function(object, value) standardGeneric("setT<-"))

## Class 'modelmoments' --------------------------------------------

#' Getter for the `mean` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getMean", function(object) standardGeneric("getMean"))

#' Getter for the `var` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getVar", function(object) standardGeneric("getVar"))

#' Getter for the `model` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getModel", function(object) standardGeneric("getModel"))

## Class 'cmodelmoments' -------------------------------------------

#' Getter for the `higher` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getHigher", function(object) standardGeneric("getHigher"))

#' Getter for the `skewness` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getSkewness", function(object) standardGeneric("getSkewness"))

#' Getter for the `kurtosis` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getKurtosis", function(object) standardGeneric("getKurtosis"))

## Class 'dmodelmoments' -------------------------------------------
#' Getter for the `over` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getOver", function(object) standardGeneric("getOver"))

#' Getter for the `factorial` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getFactorial", function(object) standardGeneric("getFactorial"))

#' Getter for `zero` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getZero", function(object) standardGeneric("getZero"))

## Class 'normultmodelmoments' -------------------------------------
#' Generates the moments of a finite mixture model
#' @export
#' @docType methods
#' @keywords internal
setGeneric("generateMoments", function(object) standardGeneric("generateMoments"))

#' Getter for the `B` slot.
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getB", function(object) standardGeneric("getB"))

#' Getter for the `W` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getW", function(object) standardGeneric("getW"))

#' Getter for the `Rdet` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getRdet", function(object) standardGeneric("getRdet"))

#' Getter for the `Rtr` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getRtr", function(object) standardGeneric("getRtr"))

#' Getter for the `corr` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getCorr", function(object) standardGeneric("getCorr"))

## Class 'exponentialmodelmoments' ---------------------------------
#' Getter for the `extrabinvar` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getExtrabinvar", function(object) standardGeneric("getExtrabinvar"))

## Class 'fdata' ----------------------------------------------------
#' Checks for the `y` slot of an `fdata` object
#' @export
#' @docType methods
#' @keywords internal
setGeneric("hasY", function(object, verbose = FALSE) standardGeneric("hasY"))

#' Checks for the `S` slot of an `fdata` object
#' @export
#' @docType methods
#' @keywords internal
setGeneric("hasS", function(object, verbose = FALSE) standardGeneric("hasS"))

#' Checks for the `exp` slot of an `fdata` object
#' @export
#' @docType methods
#' @keywords internal
setGeneric("hasExp", function(object, verbose = FALSE) standardGeneric("hasExp"))

#' Checks for the `T` slot of an `fdata` object
#' @export
#' @docType methods
#' @keywords internal
setGeneric("hasT", function(object, verbose = FALSE) standardGeneric("hasT"))

#' Getter for the `y` slot in column format
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getColY", function(object) standardGeneric("getColY"))

#' Getter for the `y` slot in row format
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getRowY", function(object) standardGeneric("getRowY"))

#' Getter for the `S` slot in column format
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getColS", function(object) standardGeneric("getColS"))

#' Getter for the `S` slot in row format
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getRowS", function(object) standardGeneric("getRowS"))

#' Getter for the `exp` slot in column format
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getColExp", function(object) standardGeneric("getColExp"))

#' Getter for the `exp` slot in row format
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getRowExp", function(object) standardGeneric("getRowExp"))

#' Getter for the `T` slot in column format
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getColT", function(object) standardGeneric("getColT"))

#' Getter for the `T` slot in row format
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getRowT", function(object) standardGeneric("getRowT"))

#' Getter for the `y` slot in stored format
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getY", function(object) standardGeneric("getY"))

#' Getter for the `bycolumn` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getBycolumn", function(object) standardGeneric("getBycolumn"))

#' Getter for the `N` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getN", function(object) standardGeneric("getN"))

#' Getter for the `S` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getS", function(object) standardGeneric("getS"))

#' Getter for the `name` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getName", function(object) standardGeneric("getName"))

#' Getter for the `type` format
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getType", function(object) standardGeneric("getType"))

#' Getter for the `sim` format
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getSim", function(object) standardGeneric("getSim"))

#' Getter for the `exp` format
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getExp", function(object) standardGeneric("getExp"))

#' Setter for the `y` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("setY<-", function(object, value) standardGeneric("setY<-"))

#' Getter for the `N` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("setN<-", function(object, value) standardGeneric("setN<-"))

#' Getter for the `S` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("setS<-", function(object, value) standardGeneric("setS<-"))

#' Setter for the `bycolumn` format
#' @export
#' @docType methods
#' @keywords internal
setGeneric("setBycolumn<-", function(object, value) standardGeneric("setBycolumn<-"))

#' Setter for the `name` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("setName<-", function(object, value) standardGeneric("setName<-"))

#' Setter for the `type` format
#' @export
#' @docType methods
#' @keywords internal
setGeneric("setType<-", function(object, value) standardGeneric("setType<-"))

#' Setter for the `sim` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("setSim<-", function(object, value) standardGeneric("setSim<-"))

#' Setter for the `exp` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("setExp<-", function(object, value) standardGeneric("setExp<-"))

## Class 'groupmoments' ----------------------------------------------
#' Getter for the `NK` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getNK", function(object) standardGeneric("getNK"))

#' Setter for the `WK` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getWK", function(object) standardGeneric("getWK"))

#' Getter for the `fdata` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getFdata", function(object) standardGeneric("getFdata"))

## Class 'sdatamoments' ----------------------------------------------
#' Getter for the `gmoments` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getGmoments", function(object) standardGeneric("getGmoments"))

## Class 'cdatamoments' ---------------------------------------------
#' Getter for the `smoments` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getSmoments", function(object) standardGeneric("getSmoments"))

## Class 'prior' -----------------------------------------------------
#' Checks for the `par` slot in the `prior` class
#' @export
#' @docType methods
#' @keywords internal
setGeneric("hasPriorPar", function(object, model, verbose = FALSE) standardGeneric("hasPriorPar"))

#' Checks for the `weight` slot in the `prior` class
#' @export
#' @docType methods
#' @keywords internal
setGeneric("hasPriorWeight", function(object, model, verbose = FALSE) standardGeneric("hasPriorWeight"))

#' Generates the prior for a specific `model` 
#' @export
#' @docType methods
#' @keywords internal
setGeneric("generatePrior", function(object, ...) standardGeneric("generatePrior"))

#' Getter for the `hier` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getHier", function(object) standardGeneric("getHier"))

#' Setter for the `hier` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("setHier<-", function(object, value) standardGeneric("setHier<-"))

## Class 'mcmc' -------------------------------------------------------
#' Getter for the `burnin` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getBurnin", function(object) standardGeneric("getBurnin"))

#' Getter for the `M` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getM", function(object) standardGeneric("getM"))

#' Getter for the `startpar` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getStartpar", function(object) standardGeneric("getStartpar"))

#' Getter for the `storeS` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getStoreS", function(object) standardGeneric("getStoreS"))

#' Getter for the `storepost` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getStorepost", function(object) standardGeneric("getStorepost"))

#' Getter for the `ranperm` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getRanperm", function(object) standardGeneric("getRanperm"))

#' Setter for the `burnin` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("setBurnin<-", function(object, value) standardGeneric("setBurnin<-"))

#' Setter for the `M` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("setM<-", function(object, value) standardGeneric("setM<-"))

#' Setter for the `startpar` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("setStartpar<-", function(object, value) standardGeneric("setStartpar<-"))

#' Setter for the `storeS` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("setStoreS<-", function(object, value) standardGeneric("setStoreS<-"))

#' Setter for the `storepost` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("setStorepost<-", function(object, value) standardGeneric("setStorepost<-"))

#' Setter for the `ranperm` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("setRanperm<-", function(object, value) standardGeneric("setRanperm<-"))

## Class 'dataclass' ----------------------------------------------------
#' Getter for the `logpy` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getLogpy", function(object) standardGeneric("getLogpy"))

#' Getter for the `prob` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getProb", function(object) standardGeneric("getProb"))

#' Getter for the mixlik slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getMixlik", function(object) standardGeneric("getMixlik"))

#' Getter for the `entropy` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getEntropy", function(object) standardGeneric("getEntropy"))

#' Getter for the `postS` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getPostS", function(object) standardGeneric("getPostS"))

#' Getter for the `loglikcd` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getLoglikcd", function(object) standardGeneric("getLoglikcd"))

## Class 'mcmcextract' --------------------------------------------------
#' Computes the model moments from MCMC samples
#' @export
#' @docType methods
#' @keywords internal
setGeneric("moments", function(object) standardGeneric("moments"))

## Class 'mcmcoutputfix' ------------------------------------------------
#' Plots the traces of the MCMC samples
#' @export
#' @docType methods
#' @keywords internal
#' @rdname plotTraces-generic
setGeneric("plotTraces", function(x, dev = TRUE, lik = 1, col = FALSE, ...) standardGeneric("plotTraces"))

#' Plots histograms of MCMC samples
#' @export
#' @docType methods
#' @keywords internal
#' @rdname plotHist-generic
setGeneric("plotHist", function(x, dev = TRUE, ...) standardGeneric("plotHist"))

#' Plots densities of MCMC samples
#' @export
#' @docType methods
#' @keywords internal
#' @rdname plotDens-generic
setGeneric("plotDens", function(x, dev = TRUE, ...) standardGeneric("plotDens"))

#' Plots sample representations of MCMC samples
#' @export
#' @docType methods
#' @keywords internal
#' @rdname plotSampRep-generic
setGeneric("plotSampRep", function(x, dev = TRUE, ...) standardGeneric("plotSampRep"))

#' Plots the posterior density of sampled component parameters
#' @export
#' @docType methods
#' @keywords internal
#' @rdname plotPostDens-generic
setGeneric("plotPostDens", function(x, dev = TRUE, ...) standardGeneric("plotPostDens"))

#' Generates a sub-chain from MCMC samples
#' @export
#' @docType methods
#' @keywords internal
#' @rdname subseq-generic
setGeneric("subseq", function(object, index) standardGeneric("subseq"))

#' Swaps elements in the MCMC sample arrays
#' @export
#' @docType methods
#' @keywords internal
#' @rdname swapElements-generic
setGeneric("swapElements", function(object, index) standardGeneric("swapElements"))

#' Extracts the MCMC samples from a specific dimension of a multivariate model
#' @export
#' @docType methods
#' @keywords internal
setGeneric("extract", function(object, index) standardGeneric("extract"))

#' Getter for the `log` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getLog", function(object) standardGeneric("getLog"))

#' Getter for the `prior` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getPrior", function(object) standardGeneric("getPrior"))

## Class 'mcmcoutputhier' -----------------------------------------------
#' Getter for the `hyper` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getHyper", function(object) standardGeneric("getHyper"))

## Class 'mcmcoutputpost' -----------------------------------------------
#' Getter for the `post` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getPost", function(object) standardGeneric("getPost"))

## Class 'mcmcoutputbase' -----------------------------------------------
#' Getter for the `ST` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getST", function(object) standardGeneric("getST"))

#' Getter for the `clust` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getClust", function(object) standardGeneric("getClust"))

## Class 'mcmcpermfix' ---------------------------------------------------
#' Getter for the `Mperm` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getMperm", function(object) standardGeneric("getMperm"))

#' Getter for the `parperm` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getParperm", function(object) standardGeneric("getParperm"))

#' Getter for the `logperm` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getLogperm", function(object) standardGeneric("getLogperm"))

## Class 'mcmcpermfixhier' -----------------------------------------------
#' Getter for the `hyperperm` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getHyperperm", function(object) standardGeneric("getHyperperm"))

## Class 'mcmcpermfixpost' -----------------------------------------------
#' Getter for the `postperm` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getPostperm", function(object) standardGeneric("getPostperm"))

## Class 'mcmcpermind' ---------------------------------------------------
#' Getter for the `relabel` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getRelabel", function(object) standardGeneric("getRelabel"))

#' Getter for the `weightperm` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getWeightperm", function(object) standardGeneric("getWeightperm"))

#' Getter for the `entropyperm` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getEntropyperm", function(object) standardGeneric("getEntropyperm"))

#' Getter for the `STperm` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getSTperm", function(object) standardGeneric("getSTperm"))

#' Getter for the `Sperm` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getSperm", function(object) standardGeneric("getSperm"))

#' Getter for the `NKperm` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getNKperm", function(object) standardGeneric("getNKperm"))

## Class 'mcmcestfix' -----------------------------------------------------
#' Getter for the `map` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getMap", function(object) standardGeneric("getMap"))

#' Getter for the `bml` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getBml", function(object) standardGeneric("getBml"))

#' Getter for the `ieavg` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getIeavg", function(object) standardGeneric("getIeavg"))

#' Getter for the `sdpost` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getSdpost", function(object) standardGeneric("getSdpost"))

## Class 'mcmcestind' ------------------------------------------------------
#' Getter for the `eavg` slot
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getEavg", function(object) standardGeneric("getEavg"))
