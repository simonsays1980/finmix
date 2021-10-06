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

## Load the dyanmic library
#' @useDynLib finmix
#' @importFrom Rcpp sourceCpp
NULL

## Class 'model' --------------------------------------------------
#' @describeIn model_class Simulates data from mixture model
setGeneric("simulate", function(model, N = 100, varargin, seed = 0) standardGeneric("simulate"))

#' @describeIn model_class Plots point process of mixture model
setGeneric("plotPointProc", function(x, dev = TRUE, ...) standardGeneric("plotPointProc"))

#' @describeIn model_class Checker for slot `weight` of model class
setGeneric("hasWeight", function(object, verbose = FALSE) standardGeneric("hasWeight"))

#' @describeIn model_class Checker for slot `T` of model class
setGeneric("hasT", function(object, verbose = FALSE) standardGeneric("hasT"))

#' @describeIn model_class Checker for slot `par` of model class
setGeneric("hasPar", function(object, verbose = FALSE) standardGeneric("hasPar"))

#' @describeIn model_class Extract marginal distribution 
setGeneric("mixturemar", function(object, J) standardGeneric("mixturemar"))

#' @describeIn model_class Getter for slot `dist` of model class
setGeneric("getDist", function(object) standardGeneric("getDist"))

#' @describeIn model_class Getter for slot `r` of model class
setGeneric("getR", function(object) standardGeneric("getR"))

#' @describeIn model_class Getter for slot `K` of model class
setGeneric("getK", function(object) standardGeneric("getK"))

#' @describeIn model_class Getter for slot `weight` of model class
setGeneric("getWeight", function(object) standardGeneric("getWeight"))

#' @describeIn model_class Getter for slot `par` of model class
setGeneric("getPar", function(object) standardGeneric("getPar"))

#' @describeIn model_class Getter for slot `indicmod` of model class
setGeneric("getIndicmod", function(object) standardGeneric("getIndicmod"))

#' @describeIn model_class Getter for slot `indicfix` of model class
setGeneric("getIndicfix", function(object) standardGeneric("getIndicfix"))

#' @describeIn model_class Getter for slot `T` of model class
setGeneric("getT", function(object) standardGeneric("getT"))

#' @describeIn model_class Setter for slot `dist` of model class
setGeneric("setDist<-", function(object, value) standardGeneric("setDist<-"))

#' @describeIn model_class Setter for slot `r` of model class
setGeneric("setR<-", function(object, value) standardGeneric("setR<-"))

#' @describeIn model_class Setter for slot `K` of model class
setGeneric("setK<-", function(object, value) standardGeneric("setK<-"))

#' @describeIn model_class Setter for slot `weight` of model class
setGeneric("setWeight<-", function(object, value) standardGeneric("setWeight<-"))

#' @describeIn model_class Setter for slot `par` of model class
setGeneric("setPar<-", function(object, value) standardGeneric("setPar<-"))

#' @describeIn model_class Setter for slot `indicmod` of model class
setGeneric("setIndicmod<-", function(object, value) standardGeneric("setIndicmod<-"))

#' @describeIn model_class Setter for slot `indicfix` of model class
setGeneric("setIndicfix<-", function(object, value) standardGeneric("setIndicfix<-"))

#' @describeIn model_class Setter for slot `T` of model class
setGeneric("setT<-", function(object, value) standardGeneric("setT<-"))

## Class 'modelmoments' --------------------------------------------

#' @describeIn modelmoments_class
setGeneric("getMean", function(object) standardGeneric("getMean"))

#' @describeIn modelmoments_class
setGeneric("getVar", function(object) standardGeneric("getVar"))

#' @describeIn modelmoments_class
setGeneric("getModel", function(object) standardGeneric("getModel"))

## Class 'cmodelmoments' -------------------------------------------

#' @describeIn modelmoments_class
setGeneric("getHigher", function(object) standardGeneric("getHigher"))

#' @describeIn modelmoments_class
setGeneric("getSkewness", function(object) standardGeneric("getSkewness"))

#' @describeIn modelmoments_class
setGeneric("getKurtosis", function(object) standardGeneric("getKurtosis"))

## Class 'dmodelmoments' -------------------------------------------
#' @describeIn modelmoments_class
setGeneric("getOver", function(object) standardGeneric("getOver"))

#' @describeIn modelmoments_class
setGeneric("getFactorial", function(object) standardGeneric("getFactorial"))

#' @describeIn modelmoments_class
setGeneric("getZero", function(object) standardGeneric("getZero"))

## Class 'normultmodelmoments' -------------------------------------
#' @describeIn modelmoments_class
setGeneric("generateMoments", function(object) standardGeneric("generateMoments"))

#' @describeIn modelmoments_class
setGeneric("getB", function(object) standardGeneric("getB"))

#' @describeIn modelmoments_class
setGeneric("getW", function(object) standardGeneric("getW"))

#' @describeIn modelmoments_class
setGeneric("getRdet", function(object) standardGeneric("getRdet"))

#' @describeIn modelmoments_class
setGeneric("getRtr", function(object) standardGeneric("getRtr"))

#' @describeIn modelmoments_class
setGeneric("getCorr", function(object) standardGeneric("getCorr"))

## Class 'exponentialmodelmoments' ---------------------------------
#' @describeIn modelmoments_class
setGeneric("getExtrabinvar", function(object) standardGeneric("getExtrabinvar"))

## Class 'fdata' ----------------------------------------------------
#' @describeIn fdata_class
setGeneric("hasY", function(object, verbose = FALSE) standardGeneric("hasY"))

#' @describeIn fdata_class
setGeneric("hasS", function(object, verbose = FALSE) standardGeneric("hasS"))

#' @describeIn fdata_class
setGeneric("hasExp", function(object, verbose = FALSE) standardGeneric("hasExp"))

#' @describeIn fdata_class
setGeneric("hasT", function(object, verbose = FALSE) standardGeneric("hasT"))

#' @describeIn fdata_class
setGeneric("getColY", function(object) standardGeneric("getColY"))

#' @describeIn fdata_class
setGeneric("getRowY", function(object) standardGeneric("getRowY"))

#' @describeIn fdata_class
setGeneric("getColS", function(object) standardGeneric("getColS"))

#' @describeIn fdata_class
setGeneric("getRowS", function(object) standardGeneric("getRowS"))

#' @describeIn fdata_class
setGeneric("getColExp", function(object) standardGeneric("getColExp"))

#' @describeIn fdata_class
setGeneric("getRowExp", function(object) standardGeneric("getRowExp"))

#' @describeIn fdata_class
setGeneric("getColT", function(object) standardGeneric("getColT"))

#' @describeIn fdata_class
setGeneric("getRowT", function(object) standardGeneric("getRowT"))

#' @describeIn fdata_class
setGeneric("getY", function(object) standardGeneric("getY"))

#' @describeIn fdata_class
setGeneric("getBycolumn", function(object) standardGeneric("getBycolumn"))

#' @describeIn fdata_class
setGeneric("getN", function(object) standardGeneric("getN"))

#' @describeIn fdata_class
setGeneric("getS", function(object) standardGeneric("getS"))

#' @describeIn fdata_class
setGeneric("getName", function(object) standardGeneric("getName"))

#' @describeIn fdata_class
setGeneric("getType", function(object) standardGeneric("getType"))

#' @describeIn fdata_class
setGeneric("getSim", function(object) standardGeneric("getSim"))

#' @describeIn fdata_class
setGeneric("getExp", function(object) standardGeneric("getExp"))

#' @describeIn fdata_class
setGeneric("setY<-", function(object, value) standardGeneric("setY<-"))

#' @describeIn fdata_class
setGeneric("setN<-", function(object, value) standardGeneric("setN<-"))

#' @describeIn fdata_class
setGeneric("setS<-", function(object, value) standardGeneric("setS<-"))

#' @describeIn fdata_class
setGeneric("setBycolumn<-", function(object, value) standardGeneric("setBycolumn<-"))

#' @describeIn fdata_class
setGeneric("setName<-", function(object, value) standardGeneric("setName<-"))

#' @describeIn fdata_class
setGeneric("setType<-", function(object, value) standardGeneric("setType<-"))

#' @describeIn fdata_class
setGeneric("setSim<-", function(object, value) standardGeneric("setSim<-"))

#' @describeIn fdata_class
setGeneric("setExp<-", function(object, value) standardGeneric("setExp<-"))

## Class 'groupmoments' ----------------------------------------------
#' @describeIn groupmoments_class
setGeneric("getNK", function(object) standardGeneric("getNK"))

#' @describeIn groupmoments_class
setGeneric("getWK", function(object) standardGeneric("getWK"))

#' @describeIn groupmoments_class
setGeneric("getFdata", function(object) standardGeneric("getFdata"))

## Class 'sdatamoments' ----------------------------------------------
#' @describeIn sdatamoments_class
setGeneric("getGmoments", function(object) standardGeneric("getGmoments"))

## Class 'cdatamoments' ---------------------------------------------
#' @describeIn cdatamoments_class
setGeneric("getSmoments", function(object) standardGeneric("getSmoments"))

## Class 'prior' -----------------------------------------------------
#' @describeIn prior-class
setGeneric("hasPriorPar", function(object, model, verbose = FALSE) standardGeneric("hasPriorPar"))

#' @describeIn prior-class
setGeneric("hasPriorWeight", function(object, model, verbose = FALSE) standardGeneric("hasPriorWeight"))

#' @describeIn prior-class
setGeneric("generatePrior", function(object, ...) standardGeneric("generatePrior"))

#' @describeIn prior-class
setGeneric("getHier", function(object) standardGeneric("getHier"))

#' @describeIn prior-class
setGeneric("setHier<-", function(object, value) standardGeneric("setHier<-"))

## Class 'mcmc' -------------------------------------------------------
#' @describeIn mcmc_class
setGeneric("getBurnin", function(object) standardGeneric("getBurnin"))

#' @describeIn mcmc_class
setGeneric("getM", function(object) standardGeneric("getM"))

#' @describeIn mcmc_class
setGeneric("getStartpar", function(object) standardGeneric("getStartpar"))

#' @describeIn mcmc_class
setGeneric("getStoreS", function(object) standardGeneric("getStoreS"))

#' @describeIn mcmc_class
setGeneric("getStorepost", function(object) standardGeneric("getStorepost"))

#' @describeIn mcmc_class
setGeneric("getRanperm", function(object) standardGeneric("getRanperm"))

#' @describeIn mcmc_class
setGeneric("setBurnin<-", function(object, value) standardGeneric("setBurnin<-"))

#' @describeIn mcmc_class
setGeneric("setM<-", function(object, value) standardGeneric("setM<-"))

#' @describeIn mcmc_class
setGeneric("setStartpar<-", function(object, value) standardGeneric("setStartpar<-"))

#' @describeIn mcmc_class
setGeneric("setStoreS<-", function(object, value) standardGeneric("setStoreS<-"))

#' @describeIn mcmc_class
setGeneric("setStorepost<-", function(object, value) standardGeneric("setStorepost<-"))

#' @describeIn mcmc_class
setGeneric("setRanperm<-", function(object, value) standardGeneric("setRanperm<-"))

## Class 'dataclass' ----------------------------------------------------
#' @describeIn dataclass
setGeneric("getLogpy", function(object) standardGeneric("getLogpy"))

#' @describeIn dataclass
setGeneric("getProb", function(object) standardGeneric("getProb"))

#' @describeIn dataclass
setGeneric("getMixlik", function(object) standardGeneric("getMixlik"))

#' @describeIn dataclass
setGeneric("getEntropy", function(object) standardGeneric("getEntropy"))

#' @describeIn dataclass
setGeneric("getPostS", function(object) standardGeneric("getPostS"))

#' @describeIn dataclass
setGeneric("getLoglikcd", function(object) standardGeneric("getLoglikcd"))

## Class 'mcmcextract' --------------------------------------------------------------------------
#' @describeIn mcmcextract_class
setGeneric("moments", function(object) standardGeneric("moments"))

## Class 'mcmcoutputfix' ------------------------------------------------
#' @describeIn mcmcoutput_class
setGeneric("plotTraces", function(x, dev = TRUE, lik = 1, col = FALSE, ...) standardGeneric("plotTraces"))

#' @describeIn mcmcoutput_class
setGeneric("plotHist", function(x, dev = TRUE, ...) standardGeneric("plotHist"))

#' @describeIn mcmcoutput_class
setGeneric("plotDens", function(x, dev = TRUE, ...) standardGeneric("plotDens"))

#' @describeIn mcmcoutput_class
setGeneric("plotSampRep", function(x, dev = TRUE, ...) standardGeneric("plotSampRep"))

#' @describeIn mcmcoutput_class
setGeneric("plotPostDens", function(x, dev = TRUE, ...) standardGeneric("plotPostDens"))

#' @describeIn mcmcoutput_class 
setGeneric("subseq", function(object, index) standardGeneric("subseq"))

#' @describeIn mcmcoutput_class
setGeneric("swapElements", function(object, index) standardGeneric("swapElements"))

#' @describeIn mcmcoutput_class
setGeneric("extract", function(object, index) standardGeneric("extract"))

#' @describeIn mcmcoutput_class
setGeneric("getLog", function(object) standardGeneric("getLog"))

#' @describeIn mcmcoutput_class
setGeneric("getPrior", function(object) standardGeneric("getPrior"))

## Class 'mcmcoutputhier' -----------------------------------------------
#' @describeIn mcmcoutput_class
setGeneric("getHyper", function(object) standardGeneric("getHyper"))

## Class 'mcmcoutputpost' -----------------------------------------------
#' @describeIn mcmcoutput_class
setGeneric("getPost", function(object) standardGeneric("getPost"))

## Class 'mcmcoutputbase' -----------------------------------------------
#' @describeIn mcmcoutput_class
setGeneric("getST", function(object) standardGeneric("getST"))

#' @describeIn mcmcoutput_class
setGeneric("getClust", function(object) standardGeneric("getClust"))

## Class 'mcmcpermfix' ---------------------------------------------------
#' @describeIn mcmcperm_class
setGeneric("getMperm", function(object) standardGeneric("getMperm"))

#' @describeIn mcmcperm_class
setGeneric("getParperm", function(object) standardGeneric("getParperm"))

#' @describeIn mcmcperm_class
setGeneric("getLogperm", function(object) standardGeneric("getLogperm"))

## Class 'mcmcpermfixhier' -----------------------------------------------
#' @noRd mcmcperm_class
setGeneric("getHyperperm", function(object) standardGeneric("getHyperperm"))

## Class 'mcmcpermfixpost' -----------------------------------------------
#' @noRd mcmcperm_class
setGeneric("getPostperm", function(object) standardGeneric("getPostperm"))

## Class 'mcmcpermind' ---------------------------------------------------
#' @describeIn mcmcperm_class
setGeneric("getRelabel", function(object) standardGeneric("getRelabel"))

#' @describeIn mcmcperm_class
setGeneric("getWeightperm", function(object) standardGeneric("getWeightperm"))

#' @describeIn mcmcperm_class
setGeneric("getEntropyperm", function(object) standardGeneric("getEntropyperm"))

#' @describeIn mcmcperm_class
setGeneric("getSTperm", function(object) standardGeneric("getSTperm"))

#' @describeIn mcmcperm_class
setGeneric("getSperm", function(object) standardGeneric("getSperm"))

#' @describeIn mcmcperm_class
setGeneric("getNKperm", function(object) standardGeneric("getNKperm"))

## Class 'mcmcestfix' -----------------------------------------------------
#' @describeIn mcmcest_class
setGeneric("getMap", function(object) standardGeneric("getMap"))

#' @describeIn mcmcest_class
setGeneric("getBml", function(object) standardGeneric("getBml"))

#' @describeIn mcmcest_class
setGeneric("getIeavg", function(object) standardGeneric("getIeavg"))

#' @describeIn mcmcest_class
setGeneric("getSdpost", function(object) standardGeneric("getSdpost"))

## Class 'mcmcestind' ------------------------------------------------------
#' @describeIn mcmcest_class
setGeneric("getEavg", function(object) standardGeneric("getEavg"))
