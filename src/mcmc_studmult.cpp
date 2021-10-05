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

//' Performs MCMC sampling for mixtures of multivariate Student-t distributions
//' 
//' @description
//' For internal usage only. This function gets passed the `fdata`, `model`,
//' `prior`, `mcmc` objects to perform MCMC sampling for a multivriate 
//' Student-t mixture model. In addition an `mcmcoutput` object is given that 
//' stores the output of MCMC sampling in R memory. Note that `finmix` relies 
//' in C++ code on so-called "mixin" layers that help to design a software by 
//' organizing code into layers that build upon each others and enable 
//' modularity in MCMC sampling by allowing to combine different sampling 
//' designs, e.g. with or without a hierarchical prior, with fixed indicators 
//' or storing posterior density parameters. See for more information on mixin 
//' layers Smaragdakis and Butory (1998). 
//' 
//' @param data_S4 An `fdata` object storing the observations and indicators.
//' @param model_S4 A `model` object specifying the multivariate Student-t 
//'   finite mixture model.
//' @param prior_S4 A `prior` object specifying the prior distribution for MCMC 
//'   sampling.
//' @param mcmc_S4 An `mcmc` object specifying the hyper-parameters for MCMC 
//'   sampling.
//' @param mcmcoutput_S4 An `mcmcoutput` object storing the outcomes from MCMC 
//'   sampling using R memory.
//' @return An `mcmcoutput` object containing the results from MCMC sampling 
//'   and using the R memory from the input argument `mcmcoutput_S4`. 
//' @export 
//' 
//' @seealso 
//' * [mixturemcmc()] for performing MCMC sampling
//' * [fdata][fdata_class] for the `fdata` class definition
//' * [model][model_class] for the `model` class definition
//' * [prior][prior_class] for the `prior` class definition
//' * [mcmc][mcmc_class] for the `mcmc` class definition
//' 
//' @references
//' * Smaragdakis, Y. and Butory, D. (1998), "Implementing layered designs with 
//'   mixin layers." In: Jul E. (eds) ECOOP’98 — Object-Oriented Programming. 
//'   ECOOP 1998. Lecture Notes in Computer Science, vol 1445. Springer, 
//'   Berlin, Heidelberg.
// [[Rcpp::export]]
RcppExport SEXP mcmc_studmult_cc(SEXP data_S4, SEXP model_S4,
                                 SEXP prior_S4, SEXP mcmc_S4, SEXP mcmcoutput_S4)
{
   /* Convert S4-classes to C++-structs */
   Rcpp::S4           dataS4O(data_S4);
   Rcpp::S4           modelS4O(model_S4);
   Rcpp::S4           priorS4O(prior_S4);
   Rcpp::S4           mcmcS4O(mcmc_S4);
   Rcpp::S4           mcmcOutputS4O(mcmcoutput_S4);
   FinmixData         finData  = FinmixData(dataS4O);
   FinmixModel        finModel = FinmixModel(modelS4O);
   FinmixPrior        finPrior = FinmixPrior(priorS4O);
   FinmixMCMC         finMCMC  = FinmixMCMC(mcmcS4O);

   const bool         INDICFIX = finModel.indicFix;
   const bool         HIER_IND = finPrior.hier;
   const bool         POST_IND = finMCMC.storePost;
   const unsigned int BURNIN   = finMCMC.burnIn;
   const unsigned int M        = finMCMC.M;
   const unsigned int K        = finModel.K;

   BASE               *ptr;

   typedef FIX<PriorStudmultFix, ParStudmultFix, LogStudmultFix,
               ParOutStudmult> STUDMULTFIX;
   typedef IND<FIX<PriorStudmultInd, ParStudmultInd, LogStudmultInd,
                   ParOutStudmult> > STUDMULTIND;
   if (INDICFIX || K == 1)
   {
      if (HIER_IND)
      {
         if (POST_IND)
         {
            ptr = new ADAPTER<POST<HIER<STUDMULTFIX,
                                        HierOutStudmult>, PostOutStudmultFix> >
                     (finData, finModel, finPrior, finMCMC,
                     mcmcOutputS4O);
         }
         else
         {
            ptr = new ADAPTER<HIER<STUDMULTFIX, HierOutStudmult> >
                     (finData, finModel, finPrior, finMCMC,
                     mcmcOutputS4O);
         }
      }
      else
      {
         if (POST_IND)
         {
            ptr = new ADAPTER<POST<STUDMULTFIX, PostOutStudmultFix> >
                     (finData, finModel, finPrior, finMCMC,
                     mcmcOutputS4O);
         }
         else
         {
            ptr = new ADAPTER<STUDMULTFIX> (finData, finModel,
                                            finPrior, finMCMC, mcmcOutputS4O);
         }
      }
   }
   else
   {
      if (HIER_IND)
      {
         if (POST_IND)
         {
            ptr = new ADAPTER<POST<HIER<STUDMULTIND,
                                        HierOutStudmult>, PostOutStudmultInd> >
                     (finData, finModel, finPrior, finMCMC,
                     mcmcOutputS4O);
         }
         else
         {
            ptr = new ADAPTER<HIER<STUDMULTIND, HierOutStudmult> >
                     (finData, finModel, finPrior, finMCMC,
                     mcmcOutputS4O);
         }
      }
      else
      {
         if (POST_IND)
         {
            ptr = new ADAPTER<POST<STUDMULTIND, PostOutStudmultInd> >
                     (finData, finModel, finPrior, finMCMC,
                     mcmcOutputS4O);
         }
         else
         {
            ptr = new ADAPTER<STUDMULTIND> (finData, finModel,
                                            finPrior, finMCMC, mcmcOutputS4O);
         }
      }
   }
   for (unsigned int i = 0; i < BURNIN + M; ++i)
   {
      ptr->update();
      ptr->store(i);
   }
   return Rcpp::wrap(mcmcOutputS4O);
}
#endif
