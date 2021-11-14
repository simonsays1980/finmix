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
* by the Free Software Foundation, either version 3 of the License, or
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

#ifndef __FINMIX_LOGBINOMIALFIX_H__
#define __FINMIX_LOGBINOMIALFIX_H__

#include "ParBinomialFix.h"
#include "likelihood.h"
#include "DataClass.h"
#include "prior_likelihood.h"

class LogBinomialFix {
public:
double mixlik;
double mixprior;

LogBinomialFix ();
virtual ~LogBinomialFix()
{
}
void update(const unsigned int&, const arma::mat&,
            const arma::ivec&, const arma::mat&, const arma::vec&,
            const ParBinomialFix&, const PriorBinomialFix&);
};

// =============================================================
// Constructor
// -------------------------------------------------------------
LogBinomialFix::LogBinomialFix () : mixlik(0.0), mixprior(0.0)
{
}

// =============================================================
// Update
// -------------------------------------------------------------

/**
 * -------------------------------------------------------------
 * update
 * @brief   Updates the log-likelihoods of the Binomial model
 * @par K           number of components
 * @par y           data matrix, N x 1
 * @par S           indicator matrix, N x 1
 * @par par         object holding the parameters
 * @par hyperPar    object holding the hyper parameters
 * @see DataClass, likelihood_binomial, priormxilik_binomial
 * @author Lars Simon Zehnder
 * -------------------------------------------------------------
 **/
void LogBinomialFix::update(const unsigned int& K, const arma::mat& y,
                            const arma::ivec& S, const arma::mat& expos, const arma::vec& T,
                            const ParBinomialFix& par, const PriorBinomialFix& hyperPar)
{
   liklist   lik   = likelihood_binomial(y, par.p, T);
   DataClass dataC = classification_fix(K, S, lik);

   mixlik = arma::sum(dataC.logLikCd);
   /* Compute likelihood of mixture prior */
   mixprior = priormixlik_binomial(par.p, hyperPar.aStart,
                                   hyperPar.bStart);
}
#endif /* __FINMIX_LOGBINOMIALFIX_H__ */



