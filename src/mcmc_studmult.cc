#ifndef __FINMIX_MCMC_STUDMULT_CC__
#define __FINMIX_MCMC_STUDMULT_CC__

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
#include "LogStudmultFix.h"
#include "ParStudmultFix.h"
#include "ParOutStudmult.h"
#include "HierOutStudmult.h"
#include "PostOutStudmultFix.h"
#include "ParStudmultInd.h"
#include "LogStudmultInd.h"
#include "PostOutStudmultInd.h"

RcppExport SEXP mcmc_studmult_cc (SEXP data_S4, SEXP model_S4,
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
    typedef FIX<PriorStudmultFix, ParStudmultFix, LogStudmultFix,
            ParOutStudmult> STUDMULTFIX;
    typedef IND<FIX<PriorStudmultInd, ParStudmultInd, LogStudmultInd,
            ParOutStudmult> > STUDMULTIND;
    if (INDICFIX || K == 1) {
        if (HIER_IND) {
            if (POST_IND) {
                ptr = new ADAPTER<POST<HIER<STUDMULTFIX,
                    HierOutStudmult>, PostOutStudmultFix> >
                        (finData, finModel, finPrior, finMCMC, 
                         mcmcOutputS4O);
            } else {
                ptr = new ADAPTER<HIER<STUDMULTFIX, HierOutStudmult> >
                    (finData, finModel, finPrior, finMCMC,
                     mcmcOutputS4O);
            }
        } else {
            if (POST_IND) {
                ptr = new ADAPTER<POST<STUDMULTFIX, PostOutStudmultFix> >
                    (finData, finModel, finPrior, finMCMC,
                     mcmcOutputS4O);
            } else {
                ptr = new ADAPTER<STUDMULTFIX> (finData, finModel,
                        finPrior, finMCMC, mcmcOutputS4O);
            }
        }
    } else {
        if (HIER_IND) {
            if (POST_IND) {                
                ptr = new ADAPTER<POST<HIER<STUDMULTIND,
                    HierOutStudmult>, PostOutStudmultInd> >
                        (finData, finModel, finPrior, finMCMC,
                         mcmcOutputS4O);
            } else {
                ptr = new ADAPTER<HIER<STUDMULTIND, HierOutStudmult> >
                    (finData, finModel, finPrior, finMCMC,
                     mcmcOutputS4O);
            }
        } else {
            if (POST_IND) {
                ptr = new ADAPTER<POST<STUDMULTIND, PostOutStudmultInd> >
                    (finData, finModel, finPrior, finMCMC,
                     mcmcOutputS4O);
            } else {
                ptr = new ADAPTER<STUDMULTIND> (finData, finModel,
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
