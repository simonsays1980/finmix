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
#ifndef __FINMIX_MCMCCONDPOISSON_CC_
#define __FINMIX_MCMCCONDPOISSON_CC_

#include <RcppArmadillo.h>              // C++ linear algebra library
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
#include "LogCondPoissonFix.h"
#include "LogCondPoissonInd.h"
#include "ParCondPoissonInd.h"
#include "ParOutCondPoisson.h"
#include "PostOutCondPoissonFix.h"
#include "PostOutCondPoissonInd.h"

//' Performs MCMC sampling for mixtures of conditional Poisson distributions
//' 
//' @description
//' For internal usage only. This function gets passed the `fdata`, `model`,
//' `prior`, `mcmc` objects to perform MCMC sampling for a conditional Poisson 
//' mixture model. In addition an `mcmcoutput` object is given that stores the 
//' output of MCMC sampling in R memory. Note that `finmix` relies in C++ code 
//' on so-called "mixin" layers that help to design a software by organizing 
//' code into layers that build upon each others and enable modularity in MCMC 
//' sampling by allowing to combine different sampling designs, e.g. with or 
//' without a hierarchical prior, with fixed indicators or storing posterior 
//' density parameters. See for more information on mixin layers Smaragdakis 
//' and Butory (1998). 
//' 
//' @param data_S4 An `fdata` object storing the observations and indicators.
//' @param model_S4 A `model` object specifying the conditional Poisson finite mixture 
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
//' * [fdata-class] for the `fdata` class definition
//' * [model-class] for the `model` class definition
//' * [prior-class] for the `prior` class definition
//' * [mcmc-class] for the `mcmc` class definition
//' 
//' @references
//' * Smaragdakis, Y. and Butory, D. (1998), "Implementing layered designs with 
//'   mixin layers." In: Jul E. (eds) ECOOP’98 — Object-Oriented Programming. 
//'   ECOOP 1998. Lecture Notes in Computer Science, vol 1445. Springer, 
//'   Berlin, Heidelberg.
// [[Rcpp::export]]
RcppExport SEXP mcmc_condpoisson_cc(SEXP data_S4, SEXP model_S4,
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
   const bool         POST_IND = finMCMC.storePost;
   const unsigned int BURNIN   = finMCMC.burnIn;
   const unsigned int M        = finMCMC.M;
   const unsigned int K        = finModel.K;

   BASE               * ptr;

   typedef FIX<PriorCondPoissonFix, ParCondPoissonFix, LogCondPoissonFix,
               ParOutCondPoisson> POISSONFIX;
   typedef IND<FIX<PriorCondPoissonInd, ParCondPoissonInd, LogCondPoissonInd,
                   ParOutCondPoisson> > POISSONIND;
   if (INDICFIX || K == 1)
   {
      if (POST_IND)
      {
         ptr = new ADAPTER<POST<POISSONFIX, PostOutCondPoissonFix> >
                  (finData, finModel, finPrior, finMCMC,
                  mcmcOutputS4O);
      }
      else
      {
         ptr = new ADAPTER<POISSONFIX> (finData, finModel,
                                        finPrior, finMCMC, mcmcOutputS4O);
      }
   }
   else
   {
      if (POST_IND)
      {
         ptr = new ADAPTER<POST<POISSONIND, PostOutCondPoissonInd> >
                  (finData, finModel, finPrior, finMCMC,
                  mcmcOutputS4O);
      }
      else
      {
         ptr = new ADAPTER<POISSONIND> (finData, finModel,
                                        finPrior, finMCMC, mcmcOutputS4O);
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
