#ifndef __FINMIX_MCMC_NORMULT_CC__
#define __FINMIX_MCMC_NORMULT_CC__

#include "FinmixData.h"
#include "FinmixModel.h"
#include "FinmixPrior.h"
#include "FinmixMCMC.h"
#include "BASE.h"
#include "ADAPTER.h"
#include "FIX.h"
#include "IND.h"
#include "HIER.h"
#include "POST.h"
#include "LogNormultFix.h"
#include "ParNormultFix.h"
#include "ParOutNormult.h"
#include "HierOutNormult.h"
#include "PostOutNormultFix.h"
#include "ParNormultInd.h"
#include "LogNormultInd.h"
#include "PostOutNormultInd.h"

//[[Rcpp::export]]
RcppExport SEXP mcmc_normult_cc (SEXP data_S4, SEXP model_S4,
        SEXP prior_S4, SEXP mcmc_S4, SEXP mcmcoutput_S4) 
{
    /* Convert S4-classes to C++-structs */
    Rcpp::S4 dataS4O(data_S4);
    Rcpp::S4 modelS4O(model_S4);
    Rcpp::S4 priorS4O(prior_S4);
    Rcpp::S4 mcmcS4O(mcmc_S4);
    Rcpp::S4 mcmcOutputS4O(mcmcoutput_S4);
    FinmixData finData      = FinmixData(dataS4O);
    FinmixModel finModel    = FinmixModel(modelS4O);
    FinmixPrior finPrior    = FinmixPrior(priorS4O);
    FinmixMCMC  finMCMC     = FinmixMCMC(mcmcS4O);

    const bool INDICFIX         = finModel.indicFix;
    const bool HIER_IND         = finPrior.hier;
    const bool POST_IND         = finMCMC.storePost;
    const unsigned int BURNIN   = finMCMC.burnIn;
    const unsigned int M        = finMCMC.M;
    const unsigned int K        = finModel.K;

    BASE *ptr;
    typedef FIX<PriorNormultFix, ParNormultFix, LogNormultFix,
            ParOutNormult> NORMULTFIX;
    typedef IND<FIX<PriorNormultInd, ParNormultInd, LogNormultInd,
            ParOutNormult> > NORMULTIND;
    if (INDICFIX || K == 1) {
        if (HIER_IND) {
            if (POST_IND) {
                ptr = new ADAPTER<POST<HIER<NORMULTFIX,
                    HierOutNormult>, PostOutNormultFix> >
                        (finData, finModel, finPrior, finMCMC, 
                         mcmcOutputS4O);
            } else {
                ptr = new ADAPTER<HIER<NORMULTFIX, HierOutNormult> >
                    (finData, finModel, finPrior, finMCMC,
                     mcmcOutputS4O);
            }
        } else {
            if (POST_IND) {
                ptr = new ADAPTER<POST<NORMULTFIX, PostOutNormultFix> >
                    (finData, finModel, finPrior, finMCMC,
                     mcmcOutputS4O);
            } else {                
                ptr = new ADAPTER<NORMULTFIX> (finData, finModel,
                        finPrior, finMCMC, mcmcOutputS4O);
            }
        }
    } else {
        if (HIER_IND) {
            if (POST_IND) {                
                ptr = new ADAPTER<POST<HIER<NORMULTIND,
                    HierOutNormult>, PostOutNormultInd> >
                        (finData, finModel, finPrior, finMCMC,
                         mcmcOutputS4O);
            } else {
                ptr = new ADAPTER<HIER<NORMULTIND, HierOutNormult> >
                    (finData, finModel, finPrior, finMCMC,
                     mcmcOutputS4O);
            }
        } else {
            if (POST_IND) {
                ptr = new ADAPTER<POST<NORMULTIND, PostOutNormultInd> >
                    (finData, finModel, finPrior, finMCMC,
                     mcmcOutputS4O);
            } else {
                ptr = new ADAPTER<NORMULTIND> (finData, finModel,
                        finPrior, finMCMC, mcmcOutputS4O);
            }
        }
    }
    for (unsigned int i = 0; i < BURNIN + M; ++i) {
        ptr->update();
        ptr->store(i);
    }
    return Rcpp::wrap(mcmcOutputS4O);
}
#endif
