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
#ifndef __FINMIX_LOGBINOMIALIND_H__
#define __FINMIX_LOGBINOMIALIND_H__

#include "LogBinomialFix.h"
#include "ParBinomialInd.h"
#include "PriorBinomialInd.h"

class LogBinomialInd : public LogBinomialFix {
public:
double cdpost;
double entropy;
double maxcdpost;

LogBinomialInd ();
virtual ~LogBinomialInd ()
{
}
void update(const unsigned int&, const arma::mat&,
            arma::ivec&, const arma::mat&, const arma::vec&,
            const ParBinomialInd&, const PriorBinomialInd&);
};

// =============================================================
// Constructor
// -------------------------------------------------------------
LogBinomialInd::LogBinomialInd () : LogBinomialFix(),
   cdpost(0.0), entropy(0.0), maxcdpost(0.0)
{
}

// =============================================================
// Update
// -------------------------------------------------------------

/**
 * -------------------------------------------------------------
 * update
 * @brief   Updates the log-likelihoods of the Binomial model
 *          and samples the indicators S.
 * @par K           number of components
 * @par y           data matrix, N x 1
 * @par S           indicator matrix from last step, N x 1
 * @par T           repetitions, N x 1
 * @par par         object holding the parameters
 * @par hyperPar    object holding the hyper parameters
 * @detail  The classification() function samples the indi-
 *          cators and computes likelihoods and entropy. As the
 *          model with fixed indicators does use a different
 *          function 'classification_fix()' it cannot be made
 *          use of inheritance, i.e. the LogBinomialFix::update()
 *          function is of no use here.
 * @see DataClass, likelihood_binomial, priormixlik_binomial,
 *      LogBinomialFix::update()
 * @author Lars Simon Zehnder
 * -------------------------------------------------------------
 **/
void LogBinomialInd::update(const unsigned int& K,
                            const arma::mat& y, arma::ivec& S, const arma::mat& expos,
                            const arma::vec& T, const ParBinomialInd& par,
                            const PriorBinomialInd& hyperPar)
{
   liklist   lik   = likelihood_binomial(y, par.p, T);
   DataClass dataC = classification(S, lik, par.weight);

   S      = dataC.newS;
   mixlik = dataC.mixLik;
   /* Compute log-likelihood of the mixture prior */
   mixprior = priormixlik_binomial(par.p, hyperPar.aStart,
                                   hyperPar.bStart);
   if (K > 1)
   {
      /* Compute log-likelihood of Dirichlet prior */
      mixprior += priormixlik_dirichlet(par.weight,
                                        hyperPar.weightStart);
      cdpost  = mixlik + mixprior + dataC.postS;
      entropy = dataC.entropy;
   }
}
#endif /* __FINMIX_LOGBINOMIALIND_H__ */



