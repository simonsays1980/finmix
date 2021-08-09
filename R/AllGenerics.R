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

## Class 'model' --------------------------------------------------
setGeneric("simulate", function(model, N = 100, varargin, seed = 0) standardGeneric("simulate"))

setGeneric("plotPointProc", function(x, dev = TRUE, ...) standardGeneric("plotPointProc"))

setGeneric("hasWeight", function(object, verbose = FALSE) standardGeneric("hasWeight"))

setGeneric("hasT", function(object, verbose = FALSE) standardGeneric("hasT"))

setGeneric("hasPar", function(object, verbose = FALSE) standardGeneric("hasPar"))

setGeneric("mixturemar", function(object, J) standardGeneric("mixturemar"))

setGeneric("getDist", function(object) standardGeneric("getDist"))

setGeneric("getR", function(object) standardGeneric("getR"))

setGeneric("getK", function(object) standardGeneric("getK"))

setGeneric("getWeight", function(object) standardGeneric("getWeight"))

setGeneric("getPar", function(object) standardGeneric("getPar"))

setGeneric("getIndicmod", function(object) standardGeneric("getIndicmod"))

setGeneric("getIndicfix", function(object) standardGeneric("getIndicfix"))

setGeneric("getT", function(object) standardGeneric("getT"))

setGeneric("setDist<-", function(object, value) standardGeneric("setDist<-"))

setGeneric("setR<-", function(object, value) standardGeneric("setR<-"))

setGeneric("setK<-", function(object, value) standardGeneric("setK<-"))

setGeneric("setWeight<-", function(object, value) standardGeneric("setWeight<-"))

setGeneric("setPar<-", function(object, value) standardGeneric("setPar<-"))

setGeneric("setIndicmod<-", function(object, value) standardGeneric("setIndicmod<-"))

setGeneric("setIndicfix<-", function(object, value) standardGeneric("setIndicfix<-"))

setGeneric("setT<-", function(object, value) standardGeneric("setT<-"))

## Class 'modelmoments' --------------------------------------------

setGeneric("getMean", function(object) standardGeneric("getMean"))

setGeneric("getVar", function(object) standardGeneric("getVar"))

setGeneric("getModel", function(object) standardGeneric("getModel"))

## Class 'cmodelmoments' -------------------------------------------

setGeneric("getHigher", function(object) standardGeneric("getHigher"))

setGeneric("getSkewness", function(object) standardGeneric("getSkewness"))

setGeneric("getKurtosis", function(object) standardGeneric("getKurtosis"))

## Class 'dmodelmoments' -------------------------------------------

setGeneric("getOver", function(object) standardGeneric("getOver"))

setGeneric("getFactorial", function(object) standardGeneric("getFactorial"))

setGeneric("getZero", function(object) standardGeneric("getZero"))

## Class 'normultmodelmoments' -------------------------------------

setGeneric("generateMoments", function(object) standardGeneric("generateMoments"))

setGeneric("getB", function(object) standardGeneric("getB"))

setGeneric("getW", function(object) standardGeneric("getW"))

setGeneric("getRdet", function(object) standardGeneric("getRdet"))

setGeneric("getRtr", function(object) standardGeneric("getRtr"))

setGeneric("getCorr", function(object) standardGeneric("getCorr"))

## Class 'exponentialmodelmoments' ---------------------------------

setGeneric("getExtrabinvar", function(object) standardGeneric("getExtrabinvar"))

## Class 'fdata' ----------------------------------------------------

setGeneric("hasY", function(object, verbose = FALSE) standardGeneric("hasY"))

setGeneric("hasS", function(object, verbose = FALSE) standardGeneric("hasS"))

setGeneric("hasExp", function(object, verbose = FALSE) standardGeneric("hasExp"))

setGeneric("hasT", function(object, verbose = FALSE) standardGeneric("hasT"))

setGeneric("getColY", function(object) standardGeneric("getColY"))

setGeneric("getRowY", function(object) standardGeneric("getRowY"))

setGeneric("getColS", function(object) standardGeneric("getColS"))

setGeneric("getRowS", function(object) standardGeneric("getRowS"))

setGeneric("getColExp", function(object) standardGeneric("getColExp"))

setGeneric("getRowExp", function(object) standardGeneric("getRowExp"))

setGeneric("getColT", function(object) standardGeneric("getColT"))

setGeneric("getRowT", function(object) standardGeneric("getRowT"))

setGeneric("getY", function(object) standardGeneric("getY"))

setGeneric("getBycolumn", function(object) standardGeneric("getBycolumn"))

setGeneric("getN", function(object) standardGeneric("getN"))

setGeneric("getS", function(object) standardGeneric("getS"))

setGeneric("getName", function(object) standardGeneric("getName"))

setGeneric("getType", function(object) standardGeneric("getType"))

setGeneric("getSim", function(object) standardGeneric("getSim"))

setGeneric("getExp", function(object) standardGeneric("getExp"))

setGeneric("setY<-", function(object, value) standardGeneric("setY<-"))

setGeneric("setN<-", function(object,  value) standardGeneric("setN<-"))

setGeneric("setS<-", function(object, value) standardGeneric("setS<-"))

setGeneric("setBycolumn<-", function(object, value) standardGeneric("setBycolumn<-"))

setGeneric("setName<-", function(object, value) standardGeneric("setName<-"))

setGeneric("setType<-", function(object, value) standardGeneric("setType<-"))

setGeneric("setSim<-", function(object, value) standardGeneric("setSim<-"))

setGeneric("setExp<-", function(object, value) standardGeneric("setExp<-"))

## Class 'groupmoments' ----------------------------------------------

setGeneric("getNK", function(object) standardGeneric("getNK"))

setGeneric("getWK", function(object) standardGeneric("getWK"))

setGeneric("getFdata", function(object) standardGeneric("getFdata"))

## Class 'sdatamoments' ----------------------------------------------

setGeneric("getGmoments", function(object) standardGeneric("getGmoments"))

## Class 'cdatamoments' ---------------------------------------------

setGeneric("getSmoments", function(object) standardGeneric("getSmoments"))

## Class 'prior' -----------------------------------------------------

setGeneric( "hasPriorPar", function( object, model, verbose = FALSE ) standardGeneric( "hasPriorPar" ) )

setGeneric("hasPriorWeight", function(object, model, verbose = FALSE) standardGeneric("hasPriorWeight"))

setGeneric("generatePrior", function(object, ... )
           {
               standardGeneric("generatePrior")
           }
)

setGeneric("getHier", function(object) standardGeneric("getHier"))

setGeneric("setHier<-", function(object, value) standardGeneric("setHier<-"))

## Class 'mcmc' -------------------------------------------------------

setGeneric("getBurnin", function(object) standardGeneric("getBurnin"))

setGeneric("getM", function(object) standardGeneric("getM"))

setGeneric("getStartpar", function(object) standardGeneric("getStartpar"))

setGeneric("getStoreS", function(object) standardGeneric("getStoreS"))

setGeneric("getStorepost", function(object) standardGeneric("getStorepost"))

setGeneric("getRanperm", function(object) standardGeneric("getRanperm"))

setGeneric("setBurnin<-", function(object, value) standardGeneric("setBurnin<-"))

setGeneric("setM<-", function(object, value) standardGeneric("setM<-"))

setGeneric("setStartpar<-", function(object, value) standardGeneric("setStartpar<-"))

setGeneric("setStoreS<-", function(object, value) standardGeneric("setStoreS<-"))

setGeneric("setStorepost<-", function(object, value) standardGeneric("setStorepost<-"))

setGeneric("setRanperm<-", function(object, value) standardGeneric("setRanperm<-"))

## Class 'dataclass' ----------------------------------------------------

setGeneric("getLogpy", function(object) standardGeneric("getLogpy"))

setGeneric("getProb", function(object) standardGeneric("getProb"))

setGeneric("getMixlik", function(object) standardGeneric("getMixlik"))

setGeneric("getEntropy", function(object) standardGeneric("getEntropy"))

setGeneric("getPostS", function(object) standardGeneric("getPostS"))

setGeneric("getLoglikcd", function(object) standardGeneric("getLoglikcd"))

## Class 'mcmcextract' --------------------------------------------------------------------------

setGeneric( "moments", function( object ) standardGeneric( "moments" ) )

## Class 'mcmcoutputfix' ------------------------------------------------

setGeneric("plotTraces", function(x, dev = TRUE, lik = 1, col = FALSE, ...) standardGeneric("plotTraces"))

setGeneric("plotHist", function(x, dev = TRUE, ...) standardGeneric("plotHist"))

setGeneric("plotDens", function(x, dev = TRUE, ...) standardGeneric("plotDens"))

setGeneric("plotSampRep", function(x, dev = TRUE, ...) standardGeneric("plotSampRep"))

setGeneric("plotPostDens", function(x, dev = TRUE, ...) standardGeneric("plotPostDens"))

setGeneric("subseq", function(object, index) standardGeneric("subseq"))

setGeneric("swapElements", function(object, index) standardGeneric("swapElements"))

setGeneric( "extract", function( object, index ) standardGeneric( "extract" ) )

setGeneric("getLog", function(object) standardGeneric("getLog"))

setGeneric("getPrior", function(object) standardGeneric("getPrior"))

## Class 'mcmcoutputhier' -----------------------------------------------

setGeneric("getHyper", function(object) standardGeneric("getHyper"))

## Class 'mcmcoutputpost' -----------------------------------------------

setGeneric("getPost", function(object) standardGeneric("getPost"))

## Class 'mcmcoutputbase' -----------------------------------------------

setGeneric("getST", function(object) standardGeneric("getST"))

setGeneric("getClust", function(object) standardGeneric("getClust"))

## Class 'mcmcpermfix' ---------------------------------------------------

setGeneric("getMperm", function(object) standardGeneric("getMperm"))

setGeneric("getParperm", function(object) standardGeneric("getParperm"))

setGeneric("getLogperm", function(object) standardGeneric("getLogperm"))

## Class 'mcmcpermfixpost' -----------------------------------------------

setGeneric("getPostperm", function(object) standardGeneric("getPostperm"))

## Class 'mcmcpermind' ---------------------------------------------------

setGeneric("getRelabel", function(object) standardGeneric("getRelabel"))

setGeneric("getWeightperm", function(object) standardGeneric("getWeightperm"))

setGeneric("getEntropyperm", function(object) standardGeneric("getEntropyperm"))

setGeneric("getSTperm", function(object) standardGeneric("getSTperm"))

setGeneric("getSperm", function(object) standardGeneric("getSperm"))

setGeneric("getNKperm", function(object) standardGeneric("getNKperm"))

## Class 'mcmcestfix' -----------------------------------------------------

setGeneric("getMap", function(object) standardGeneric("getMap"))

setGeneric("getBml", function(object) standardGeneric("getBml"))

setGeneric("getIeavg", function(object) standardGeneric("getIeavg"))

setGeneric("getSdpost", function(object) standardGeneric("getSdpost"))

## Class 'mcmcestind' ------------------------------------------------------

setGeneric("getEavg", function(object) standardGeneric("getEavg"))













