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
#' @noRd
setGeneric("simulate", function(model, N = 100, varargin, seed = 0) standardGeneric("simulate"))

#' @noRd
setGeneric("plotPointProc", function(x, dev = TRUE, ...) standardGeneric("plotPointProc"))

#' @noRd
setGeneric("hasWeight", function(object, verbose = FALSE) standardGeneric("hasWeight"))

#' @noRd
setGeneric("hasT", function(object, verbose = FALSE) standardGeneric("hasT"))

#' @noRd
setGeneric("hasPar", function(object, verbose = FALSE) standardGeneric("hasPar"))

#' @noRd 
setGeneric("mixturemar", function(object, J) standardGeneric("mixturemar"))

#' @noRd
setGeneric("getDist", function(object) standardGeneric("getDist"))

#' @noRd
setGeneric("getR", function(object) standardGeneric("getR"))

#' @noRd
setGeneric("getK", function(object) standardGeneric("getK"))

#' @noRd
setGeneric("getWeight", function(object) standardGeneric("getWeight"))

#' @noRd
setGeneric("getPar", function(object) standardGeneric("getPar"))

#' @noRd
setGeneric("getIndicmod", function(object) standardGeneric("getIndicmod"))

#' @noRd
setGeneric("getIndicfix", function(object) standardGeneric("getIndicfix"))

#' @noRd
setGeneric("getT", function(object) standardGeneric("getT"))

#' @noRd
setGeneric("setDist<-", function(object, value) standardGeneric("setDist<-"))

#' @noRd
setGeneric("setR<-", function(object, value) standardGeneric("setR<-"))

#' @noRd
setGeneric("setK<-", function(object, value) standardGeneric("setK<-"))

#' @noRd
setGeneric("setWeight<-", function(object, value) standardGeneric("setWeight<-"))

#' @noRd
setGeneric("setPar<-", function(object, value) standardGeneric("setPar<-"))

#' @noRd
setGeneric("setIndicmod<-", function(object, value) standardGeneric("setIndicmod<-"))

#' @noRd
setGeneric("setIndicfix<-", function(object, value) standardGeneric("setIndicfix<-"))

#' @noRd
setGeneric("setT<-", function(object, value) standardGeneric("setT<-"))

## Class 'modelmoments' --------------------------------------------

#' @noRd
setGeneric("getMean", function(object) standardGeneric("getMean"))

#' @noRd
setGeneric("getVar", function(object) standardGeneric("getVar"))

#' @noRd
setGeneric("getModel", function(object) standardGeneric("getModel"))

## Class 'cmodelmoments' -------------------------------------------

#' @noRd
setGeneric("getHigher", function(object) standardGeneric("getHigher"))

#' @noRd
setGeneric("getSkewness", function(object) standardGeneric("getSkewness"))

#' @noRd
setGeneric("getKurtosis", function(object) standardGeneric("getKurtosis"))

## Class 'dmodelmoments' -------------------------------------------
#' @noRd
setGeneric("getOver", function(object) standardGeneric("getOver"))

#' @noRd
setGeneric("getFactorial", function(object) standardGeneric("getFactorial"))

#' @noRd
setGeneric("getZero", function(object) standardGeneric("getZero"))

## Class 'normultmodelmoments' -------------------------------------
#' @noRd
setGeneric("generateMoments", function(object) standardGeneric("generateMoments"))

#' @noRd
setGeneric("getB", function(object) standardGeneric("getB"))

#' @noRd
setGeneric("getW", function(object) standardGeneric("getW"))

#' @noRd
setGeneric("getRdet", function(object) standardGeneric("getRdet"))

#' @noRd
setGeneric("getRtr", function(object) standardGeneric("getRtr"))

#' @noRd
setGeneric("getCorr", function(object) standardGeneric("getCorr"))

## Class 'exponentialmodelmoments' ---------------------------------
#' @noRd
setGeneric("getExtrabinvar", function(object) standardGeneric("getExtrabinvar"))

## Class 'fdata' ----------------------------------------------------
#' @noRd
setGeneric("hasY", function(object, verbose = FALSE) standardGeneric("hasY"))

#' @noRd
setGeneric("hasS", function(object, verbose = FALSE) standardGeneric("hasS"))

#' @noRd
setGeneric("hasExp", function(object, verbose = FALSE) standardGeneric("hasExp"))

#' @noRd
setGeneric("hasT", function(object, verbose = FALSE) standardGeneric("hasT"))

#' @noRd
setGeneric("getColY", function(object) standardGeneric("getColY"))

#' @noRd
setGeneric("getRowY", function(object) standardGeneric("getRowY"))

#' @noRd
setGeneric("getColS", function(object) standardGeneric("getColS"))

#' @noRd
setGeneric("getRowS", function(object) standardGeneric("getRowS"))

#' @noRd
setGeneric("getColExp", function(object) standardGeneric("getColExp"))

#' @noRd
setGeneric("getRowExp", function(object) standardGeneric("getRowExp"))

#' @noRd
setGeneric("getColT", function(object) standardGeneric("getColT"))

#' @noRd
setGeneric("getRowT", function(object) standardGeneric("getRowT"))

#' @noRd
setGeneric("getY", function(object) standardGeneric("getY"))

#' @noRd
setGeneric("getBycolumn", function(object) standardGeneric("getBycolumn"))

#' @noRd
setGeneric("getN", function(object) standardGeneric("getN"))

#' @noRd
setGeneric("getS", function(object) standardGeneric("getS"))

#' @noRd
setGeneric("getName", function(object) standardGeneric("getName"))

#' @noRd
setGeneric("getType", function(object) standardGeneric("getType"))

#' @noRd
setGeneric("getSim", function(object) standardGeneric("getSim"))

#' @noRd
setGeneric("getExp", function(object) standardGeneric("getExp"))

#' @noRd
setGeneric("setY<-", function(object, value) standardGeneric("setY<-"))

#' @noRd
setGeneric("setN<-", function(object, value) standardGeneric("setN<-"))

#' @noRd
setGeneric("setS<-", function(object, value) standardGeneric("setS<-"))

#' @noRd
setGeneric("setBycolumn<-", function(object, value) standardGeneric("setBycolumn<-"))

#' @noRd
setGeneric("setName<-", function(object, value) standardGeneric("setName<-"))

#' @noRd
setGeneric("setType<-", function(object, value) standardGeneric("setType<-"))

#' @noRd
setGeneric("setSim<-", function(object, value) standardGeneric("setSim<-"))

#' @noRd
setGeneric("setExp<-", function(object, value) standardGeneric("setExp<-"))

## Class 'groupmoments' ----------------------------------------------
#' @noRd
setGeneric("getNK", function(object) standardGeneric("getNK"))

#' @noRd
setGeneric("getWK", function(object) standardGeneric("getWK"))

#' @noRd
setGeneric("getFdata", function(object) standardGeneric("getFdata"))

## Class 'sdatamoments' ----------------------------------------------
#' @noRd
setGeneric("getGmoments", function(object) standardGeneric("getGmoments"))

## Class 'cdatamoments' ---------------------------------------------
#' @noRd
setGeneric("getSmoments", function(object) standardGeneric("getSmoments"))

## Class 'prior' -----------------------------------------------------
#' @noRd
setGeneric("hasPriorPar", function(object, model, verbose = FALSE) standardGeneric("hasPriorPar"))

#' @noRd
setGeneric("hasPriorWeight", function(object, model, verbose = FALSE) standardGeneric("hasPriorWeight"))

#' @noRd
setGeneric("generatePrior", function(object, ...) standardGeneric("generatePrior"))

#' @noRd
setGeneric("getHier", function(object) standardGeneric("getHier"))

#' @noRd
setGeneric("setHier<-", function(object, value) standardGeneric("setHier<-"))

## Class 'mcmc' -------------------------------------------------------
#' @noRd
setGeneric("getBurnin", function(object) standardGeneric("getBurnin"))

#' @noRd
setGeneric("getM", function(object) standardGeneric("getM"))

#' @noRd
setGeneric("getStartpar", function(object) standardGeneric("getStartpar"))

#' @noRd
setGeneric("getStoreS", function(object) standardGeneric("getStoreS"))

#' @noRd
setGeneric("getStorepost", function(object) standardGeneric("getStorepost"))

#' @noRd
setGeneric("getRanperm", function(object) standardGeneric("getRanperm"))

#' @noRd
setGeneric("setBurnin<-", function(object, value) standardGeneric("setBurnin<-"))

#' @noRd
setGeneric("setM<-", function(object, value) standardGeneric("setM<-"))

#' @noRd
setGeneric("setStartpar<-", function(object, value) standardGeneric("setStartpar<-"))

#' @noRd
setGeneric("setStoreS<-", function(object, value) standardGeneric("setStoreS<-"))

#' @noRd
setGeneric("setStorepost<-", function(object, value) standardGeneric("setStorepost<-"))

#' @noRd
setGeneric("setRanperm<-", function(object, value) standardGeneric("setRanperm<-"))

## Class 'dataclass' ----------------------------------------------------
#' @noRd
setGeneric("getLogpy", function(object) standardGeneric("getLogpy"))

#' @noRd
setGeneric("getProb", function(object) standardGeneric("getProb"))

#' @noRd
setGeneric("getMixlik", function(object) standardGeneric("getMixlik"))

#' @noRd
setGeneric("getEntropy", function(object) standardGeneric("getEntropy"))

#' @noRd
setGeneric("getPostS", function(object) standardGeneric("getPostS"))

#' @noRd
setGeneric("getLoglikcd", function(object) standardGeneric("getLoglikcd"))

## Class 'mcmcextract' --------------------------------------------------------------------------
#' @noRd
setGeneric("moments", function(object) standardGeneric("moments"))

## Class 'mcmcoutputfix' ------------------------------------------------
#' @noRd
setGeneric("plotTraces", function(x, dev = TRUE, lik = 1, col = FALSE, ...) standardGeneric("plotTraces"))

#' @noRd
setGeneric("plotHist", function(x, dev = TRUE, ...) standardGeneric("plotHist"))

#' @noRd
setGeneric("plotDens", function(x, dev = TRUE, ...) standardGeneric("plotDens"))

#' @noRd
setGeneric("plotSampRep", function(x, dev = TRUE, ...) standardGeneric("plotSampRep"))

#' @noRd
setGeneric("plotPostDens", function(x, dev = TRUE, ...) standardGeneric("plotPostDens"))

#' @noRd
setGeneric("subseq", function(object, index) standardGeneric("subseq"))

#' @noRd
setGeneric("swapElements", function(object, index) standardGeneric("swapElements"))

#' @noRd
setGeneric("extract", function(object, index) standardGeneric("extract"))

#' @noRd
setGeneric("getLog", function(object) standardGeneric("getLog"))

#' @noRd
setGeneric("getPrior", function(object) standardGeneric("getPrior"))

## Class 'mcmcoutputhier' -----------------------------------------------
#' @noRd
setGeneric("getHyper", function(object) standardGeneric("getHyper"))

## Class 'mcmcoutputpost' -----------------------------------------------
#' @noRd
setGeneric("getPost", function(object) standardGeneric("getPost"))

## Class 'mcmcoutputbase' -----------------------------------------------
#' @noRd
setGeneric("getST", function(object) standardGeneric("getST"))

#' @noRd
setGeneric("getClust", function(object) standardGeneric("getClust"))

## Class 'mcmcpermfix' ---------------------------------------------------
#' @noRd
setGeneric("getMperm", function(object) standardGeneric("getMperm"))

#' @noRd
setGeneric("getParperm", function(object) standardGeneric("getParperm"))

#' @noRd
setGeneric("getLogperm", function(object) standardGeneric("getLogperm"))

## Class 'mcmcpermfixhier' -----------------------------------------------
#' @noRd
setGeneric("getHyperperm", function(object) standardGeneric("getHyperperm"))

## Class 'mcmcpermfixpost' -----------------------------------------------
#' @noRd
setGeneric("getPostperm", function(object) standardGeneric("getPostperm"))

## Class 'mcmcpermind' ---------------------------------------------------
#' @noRd
setGeneric("getRelabel", function(object) standardGeneric("getRelabel"))

#' @noRd
setGeneric("getWeightperm", function(object) standardGeneric("getWeightperm"))

#' @noRd
setGeneric("getEntropyperm", function(object) standardGeneric("getEntropyperm"))

#' @noRd
setGeneric("getSTperm", function(object) standardGeneric("getSTperm"))

#' @noRd
setGeneric("getSperm", function(object) standardGeneric("getSperm"))

#' @noRd
setGeneric("getNKperm", function(object) standardGeneric("getNKperm"))

## Class 'mcmcestfix' -----------------------------------------------------
#' @noRd
setGeneric("getMap", function(object) standardGeneric("getMap"))

#' @noRd
setGeneric("getBml", function(object) standardGeneric("getBml"))

#' @noRd
setGeneric("getIeavg", function(object) standardGeneric("getIeavg"))

#' @noRd
setGeneric("getSdpost", function(object) standardGeneric("getSdpost"))

## Class 'mcmcestind' ------------------------------------------------------
#' @noRd
setGeneric("getEavg", function(object) standardGeneric("getEavg"))
