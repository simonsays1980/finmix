/******************************************************************************
*
* Copyright (C) 2013 Lars Simon Zehnder. All Rights Reserved.
*
* Author: Lars Simon Zehnder <simon.zehnder@gmail.com>
*
* This file is part of the R package finmix.
*
* finmix is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published
* by the Free Software Foundatio, either version 3 of the License, or
* any later version.
*
* finmix is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with finmix. If not, see <http://www.gnu.org/licenses/>.
*
******************************************************************************/
#include "FinmixData.h"
#include "FinmixModel.h"
#include "FinmixPrior.h"
#include "FinmixMCMC.h"
#include "ADAPTER.h"
#include "FIX.h"
#include "IND.h"
#include "POST.h"
#include "PriorBinomialFix.h"
#include "ParBinomialFix.h"
#include "PriorBinomialInd.h"
#include "ParBinomialInd.h"
#include "LogBinomialInd.h"
#include "ParOutBinomial.h"
#include "PostOutBinomialInd.h"

//' Performs MCMC sampling for mixtures of Binomial distributions
//' 
//' @description
//' For internal usage only. This function gets passed the `fdata`, `model`,
//' `prior`, `mcmc` objects to perform MCMC sampling for a Binomial mixture 
//' model. In addition an `mcmcoutput` object is given that stores the output 
//' of MCMC sampling in R memory. Note that `finmix` relies in C++ code on 
//' so-called "mixin" layers that help to design a software by organizing code 
//' into layers that build upon each others and enable modularity in MCMC 
//' sampling by allowing to combine different sampling designs, e.g. with or 
//' without a hierarchical prior, with fixed indicators or storing posterior 
//' density parameters. See for more information on mixin layers Smaragdakis 
//' and Butory (1998). 
//' 
//' @param data_S4 An `fdata` object storing the observations and indicators.
//' @param model_S4 A `model` object specifying the Binomial finite mixture 
//'   model.
//' @param prior_S4 A `prior` object specifying the prior distribution for MCMC 
//'   sampling.
//' @param mcmc_S4 A `mcmc` object specifying the hyper-parameters for MCMC 
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
RcppExport SEXP mcmc_binomial_cc(SEXP fdata_S4, SEXP model_S4,
                                 SEXP prior_S4, SEXP mcmc_S4, SEXP mcmcoutput_S4)
{
   /* Convert S4-classes to C++-classes */
   Rcpp::S4           fdataS4(fdata_S4);
   Rcpp::S4           modelS4(model_S4);
   Rcpp::S4           priorS4(prior_S4);
   Rcpp::S4           mcmcS4(mcmc_S4);
   Rcpp::S4           mcmcoutputS4(mcmcoutput_S4);
   FinmixData         finFdata = FinmixData(fdataS4);
   FinmixModel        finModel = FinmixModel(modelS4);
   FinmixPrior        finPrior = FinmixPrior(priorS4);
   FinmixMCMC         finMCMC  = FinmixMCMC(mcmcS4);

   const bool         INDICFIX = finModel.indicFix;
   const bool         POST_IND = finMCMC.storePost;
   const unsigned int BURNIN   = finMCMC.burnIn;
   const unsigned int M        = finMCMC.M;
   const unsigned int K        = finModel.K;

   BASE               * ptr;

   typedef FIX<PriorBinomialFix, ParBinomialFix, LogBinomialFix,
               ParOutBinomial> BINOMIALFIX;
   typedef IND<FIX<PriorBinomialInd, ParBinomialInd, LogBinomialInd,
                   ParOutBinomial> > BINOMIALIND;
   if (INDICFIX || K == 1)
   {
      if (POST_IND)
      {
         ptr = new ADAPTER<POST<BINOMIALFIX, PostOutBinomialFix> >
                  (finFdata, finModel, finPrior, finMCMC,
                  mcmcoutputS4);
      }
      else
      {
         ptr = new ADAPTER<BINOMIALFIX>
                  (finFdata, finModel, finPrior, finMCMC,
                  mcmcoutputS4);
      }
   }
   else
   {
      if (POST_IND)
      {
         ptr = new ADAPTER<POST<BINOMIALIND, PostOutBinomialInd> >
                  (finFdata, finModel, finPrior, finMCMC,
                  mcmcoutputS4);
      }
      else
      {
         ptr = new ADAPTER<BINOMIALIND>
                  (finFdata, finModel, finPrior, finMCMC,
                  mcmcoutputS4);
      }
   }
   for (unsigned int i = 0; i < BURNIN + M; ++i)
   {
      ptr->update();
      ptr->store(i);
   }
   return Rcpp::wrap(mcmcoutputS4);
}
