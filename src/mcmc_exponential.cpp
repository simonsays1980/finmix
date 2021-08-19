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
#ifndef _FINMIX_MCMCEXPONENTIAL_CC__
#define _FINMIX_MCMCEXPONENTIAL_CC__

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
#include "LogExponentialFix.h"
#include "LogExponentialInd.h"
#include "ParExponentialInd.h"
#include "ParOutExponential.h"
#include "PostOutExponentialFix.h"
#include "PostOutExponentialInd.h"


RcppExport SEXP mcmc_exponential_cc(SEXP data_S4, SEXP model_S4, 
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
   	FinmixMCMC finMCMC      = FinmixMCMC(mcmcS4O);

	const bool INDICFIX         = finModel.indicFix;
	const bool HIER_IND         = finPrior.hier;
	const bool POST_IND         = finMCMC.storePost;
	const unsigned int BURNIN   = finMCMC.burnIn;
	const unsigned int M 		= finMCMC.M;
	const unsigned int K 		= finModel.K;
	
	BASE* ptr;
	typedef FIX<PriorExponentialFix, ParExponentialFix, LogExponentialFix, 
		ParOutExponential> EXPONENTIALFIX;
	typedef IND<FIX<PriorExponentialInd, ParExponentialInd, LogExponentialInd, 
		ParOutExponential> > EXPONENTIALIND;
	if (INDICFIX || K == 1) 
	{ 
        if (POST_IND) {
            ptr = new ADAPTER<POST<EXPONENTIALFIX, PostOutExponentialFix> >
                (finData, finModel, finPrior, finMCMC, 
                mcmcOutputS4O);
        }
        else {
            ptr = new ADAPTER<EXPONENTIALFIX> (finData, finModel,
                finPrior, finMCMC, mcmcOutputS4O);
        }
	}
	else {	
        if (POST_IND) {
            ptr = new ADAPTER<POST<EXPONENTIALIND, PostOutExponentialInd> >
                (finData, finModel, finPrior, finMCMC,
                mcmcOutputS4O);
        }
        else {
            ptr = new ADAPTER<EXPONENTIALIND> (finData, finModel,
                finPrior, finMCMC, mcmcOutputS4O);
        }
    }
	
	for(unsigned int i = 0; i < BURNIN + M; ++i) {
		ptr->update();
		ptr->store(i);
	}		
		
	return Rcpp::wrap(mcmcOutputS4O);	
}
#endif // __FINMIX_MCMCEXPONENTIAL_CC__
