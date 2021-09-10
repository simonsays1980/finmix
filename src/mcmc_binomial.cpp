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
