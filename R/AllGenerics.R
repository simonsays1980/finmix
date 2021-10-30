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
#' @export
#' @docType methods
#' @keywords internal
setGeneric("simulate", function(model, N = 100, varargin, seed = 0) standardGeneric("simulate"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("plotPointProc", function(x, dev = TRUE, ...) standardGeneric("plotPointProc"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("hasWeight", function(object, verbose = FALSE) standardGeneric("hasWeight"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("hasT", function(object, verbose = FALSE) standardGeneric("hasT"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("hasPar", function(object, verbose = FALSE) standardGeneric("hasPar"))

#' @export
#' @docType methods
#' @keywords internal 
setGeneric("mixturemar", function(object, J) standardGeneric("mixturemar"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getDist", function(object) standardGeneric("getDist"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getR", function(object) standardGeneric("getR"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getK", function(object) standardGeneric("getK"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getWeight", function(object) standardGeneric("getWeight"))

#' @export 
#' @docType methods
#' @keywords internal
setGeneric("getPar", function(object) standardGeneric("getPar"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getIndicmod", function(object) standardGeneric("getIndicmod"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getIndicfix", function(object) standardGeneric("getIndicfix"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getT", function(object) standardGeneric("getT"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("setDist<-", function(object, value) standardGeneric("setDist<-"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("setR<-", function(object, value) standardGeneric("setR<-"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("setK<-", function(object, value) standardGeneric("setK<-"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("setWeight<-", function(object, value) standardGeneric("setWeight<-"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("setPar<-", function(object, value) standardGeneric("setPar<-"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("setIndicmod<-", function(object, value) standardGeneric("setIndicmod<-"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("setIndicfix<-", function(object, value) standardGeneric("setIndicfix<-"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("setT<-", function(object, value) standardGeneric("setT<-"))

## Class 'modelmoments' --------------------------------------------

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getMean", function(object) standardGeneric("getMean"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getVar", function(object) standardGeneric("getVar"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getModel", function(object) standardGeneric("getModel"))

## Class 'cmodelmoments' -------------------------------------------

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getHigher", function(object) standardGeneric("getHigher"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getSkewness", function(object) standardGeneric("getSkewness"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getKurtosis", function(object) standardGeneric("getKurtosis"))

## Class 'dmodelmoments' -------------------------------------------
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getOver", function(object) standardGeneric("getOver"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getFactorial", function(object) standardGeneric("getFactorial"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getZero", function(object) standardGeneric("getZero"))

## Class 'normultmodelmoments' -------------------------------------
#' @export
#' @docType methods
#' @keywords internal
#' @aliases generateMoments,cmodelmoments-class,csmodelmoments-class,
#' exponentialmodelmoments,binomialmodelmoments-class
setGeneric("generateMoments", function(object) standardGeneric("generateMoments"))

#' @export
#' @docType methods
#' @keywords internal
#' @aliases getB,cmodelmoments-class,exponentialmodelmoments-class
setGeneric("getB", function(object) standardGeneric("getB"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getW", function(object) standardGeneric("getW"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getRdet", function(object) standardGeneric("getRdet"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getRtr", function(object) standardGeneric("getRtr"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getCorr", function(object) standardGeneric("getCorr"))

## Class 'exponentialmodelmoments' ---------------------------------
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getExtrabinvar", function(object) standardGeneric("getExtrabinvar"))

## Class 'fdata' ----------------------------------------------------
#' @export
#' @docType methods
#' @keywords internal
setGeneric("hasY", function(object, verbose = FALSE) standardGeneric("hasY"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("hasS", function(object, verbose = FALSE) standardGeneric("hasS"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("hasExp", function(object, verbose = FALSE) standardGeneric("hasExp"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("hasT", function(object, verbose = FALSE) standardGeneric("hasT"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getColY", function(object) standardGeneric("getColY"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getRowY", function(object) standardGeneric("getRowY"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getColS", function(object) standardGeneric("getColS"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getRowS", function(object) standardGeneric("getRowS"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getColExp", function(object) standardGeneric("getColExp"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getRowExp", function(object) standardGeneric("getRowExp"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getColT", function(object) standardGeneric("getColT"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getRowT", function(object) standardGeneric("getRowT"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getY", function(object) standardGeneric("getY"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getBycolumn", function(object) standardGeneric("getBycolumn"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getN", function(object) standardGeneric("getN"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getS", function(object) standardGeneric("getS"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getName", function(object) standardGeneric("getName"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getType", function(object) standardGeneric("getType"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getSim", function(object) standardGeneric("getSim"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getExp", function(object) standardGeneric("getExp"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("setY<-", function(object, value) standardGeneric("setY<-"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("setN<-", function(object, value) standardGeneric("setN<-"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("setS<-", function(object, value) standardGeneric("setS<-"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("setBycolumn<-", function(object, value) standardGeneric("setBycolumn<-"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("setName<-", function(object, value) standardGeneric("setName<-"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("setType<-", function(object, value) standardGeneric("setType<-"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("setSim<-", function(object, value) standardGeneric("setSim<-"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("setExp<-", function(object, value) standardGeneric("setExp<-"))

## Class 'groupmoments' ----------------------------------------------
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getNK", function(object) standardGeneric("getNK"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getWK", function(object) standardGeneric("getWK"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getFdata", function(object) standardGeneric("getFdata"))

## Class 'sdatamoments' ----------------------------------------------
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getGmoments", function(object) standardGeneric("getGmoments"))

## Class 'cdatamoments' ---------------------------------------------
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getSmoments", function(object) standardGeneric("getSmoments"))

## Class 'prior' -----------------------------------------------------
#' @export
#' @docType methods
#' @keywords internal
setGeneric("hasPriorPar", function(object, model, verbose = FALSE) standardGeneric("hasPriorPar"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("hasPriorWeight", function(object, model, verbose = FALSE) standardGeneric("hasPriorWeight"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("generatePrior", function(object, ...) standardGeneric("generatePrior"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getHier", function(object) standardGeneric("getHier"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("setHier<-", function(object, value) standardGeneric("setHier<-"))

## Class 'mcmc' -------------------------------------------------------
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getBurnin", function(object) standardGeneric("getBurnin"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getM", function(object) standardGeneric("getM"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getStartpar", function(object) standardGeneric("getStartpar"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getStoreS", function(object) standardGeneric("getStoreS"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getStorepost", function(object) standardGeneric("getStorepost"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getRanperm", function(object) standardGeneric("getRanperm"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("setBurnin<-", function(object, value) standardGeneric("setBurnin<-"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("setM<-", function(object, value) standardGeneric("setM<-"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("setStartpar<-", function(object, value) standardGeneric("setStartpar<-"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("setStoreS<-", function(object, value) standardGeneric("setStoreS<-"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("setStorepost<-", function(object, value) standardGeneric("setStorepost<-"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("setRanperm<-", function(object, value) standardGeneric("setRanperm<-"))

## Class 'dataclass' ----------------------------------------------------
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getLogpy", function(object) standardGeneric("getLogpy"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getProb", function(object) standardGeneric("getProb"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getMixlik", function(object) standardGeneric("getMixlik"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getEntropy", function(object) standardGeneric("getEntropy"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getPostS", function(object) standardGeneric("getPostS"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getLoglikcd", function(object) standardGeneric("getLoglikcd"))

## Class 'mcmcextract' --------------------------------------------------------------------------
#' @export
#' @docType methods
#' @keywords internal
setGeneric("moments", function(object) standardGeneric("moments"))

## Class 'mcmcoutputfix' ------------------------------------------------
#' @export
#' @docType methods
#' @keywords internal
setGeneric("plotTraces", function(x, dev = TRUE, lik = 1, col = FALSE, ...) standardGeneric("plotTraces"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("plotHist", function(x, dev = TRUE, ...) standardGeneric("plotHist"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("plotDens", function(x, dev = TRUE, ...) standardGeneric("plotDens"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("plotSampRep", function(x, dev = TRUE, ...) standardGeneric("plotSampRep"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("plotPostDens", function(x, dev = TRUE, ...) standardGeneric("plotPostDens"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("subseq", function(object, index) standardGeneric("subseq"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("swapElements", function(object, index) standardGeneric("swapElements"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("extract", function(object, index) standardGeneric("extract"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getLog", function(object) standardGeneric("getLog"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getPrior", function(object) standardGeneric("getPrior"))

## Class 'mcmcoutputhier' -----------------------------------------------
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getHyper", function(object) standardGeneric("getHyper"))

## Class 'mcmcoutputpost' -----------------------------------------------
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getPost", function(object) standardGeneric("getPost"))

## Class 'mcmcoutputbase' -----------------------------------------------
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getST", function(object) standardGeneric("getST"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getClust", function(object) standardGeneric("getClust"))

## Class 'mcmcpermfix' ---------------------------------------------------
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getMperm", function(object) standardGeneric("getMperm"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getParperm", function(object) standardGeneric("getParperm"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getLogperm", function(object) standardGeneric("getLogperm"))

## Class 'mcmcpermfixhier' -----------------------------------------------
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getHyperperm", function(object) standardGeneric("getHyperperm"))

## Class 'mcmcpermfixpost' -----------------------------------------------
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getPostperm", function(object) standardGeneric("getPostperm"))

## Class 'mcmcpermind' ---------------------------------------------------
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getRelabel", function(object) standardGeneric("getRelabel"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getWeightperm", function(object) standardGeneric("getWeightperm"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getEntropyperm", function(object) standardGeneric("getEntropyperm"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getSTperm", function(object) standardGeneric("getSTperm"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getSperm", function(object) standardGeneric("getSperm"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getNKperm", function(object) standardGeneric("getNKperm"))

## Class 'mcmcestfix' -----------------------------------------------------
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getMap", function(object) standardGeneric("getMap"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getBml", function(object) standardGeneric("getBml"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getIeavg", function(object) standardGeneric("getIeavg"))

#' @export
#' @docType methods
#' @keywords internal
setGeneric("getSdpost", function(object) standardGeneric("getSdpost"))

## Class 'mcmcestind' ------------------------------------------------------
#' @export
#' @docType methods
#' @keywords internal
setGeneric("getEavg", function(object) standardGeneric("getEavg"))
