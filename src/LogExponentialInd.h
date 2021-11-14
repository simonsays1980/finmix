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
* along with 'finmix'. If not, see <http://www.gnu.org/licenses/>.
*
******************************************************************************/
#ifndef __FINMIX_LOGEXPONENTIALIND_H__
#define __FINMIX_LOGEXPONENTIALIND_H__

#include "LogExponentialFix.h"
#include "ParExponentialInd.h"
#include "PriorExponentialInd.h"

class LogExponentialInd : public LogExponentialFix {
public:
double cdpost;
double entropy;
double maxcdpost;

LogExponentialInd ();
virtual ~LogExponentialInd ()
{
}
void update(const unsigned int&, const arma::mat&,
            arma::ivec&, const arma::mat&, const arma::vec&,
            const ParExponentialInd&, const PriorExponentialInd&);
};

// =============================================================
// Constructor
// -------------------------------------------------------------
LogExponentialInd::LogExponentialInd () : LogExponentialFix(),
   cdpost(0.0), entropy(0.0), maxcdpost(0.0)
{
}

/**
 * -------------------------------------------------------------
 * update
 * @brief   Updates the log-likelihoods of the Exponential model
 *          and samples the indicators S.
 * @par K           number of components
 * @par y           data matrix, N x 1
 * @par S           indicator matrix from last step, N x 1
 * @par par         object holding the parameters
 * @par hyperPar    object holding the hyper parameters
 * @detail  The classification() function samples the indi-
 *          cators and computes likelihoods and entropy. As the
 *          model with fixed indicators does use a different
 *          function 'classification_fix()' it cannot be made
 *          use of inheritance, i.e. the LogExponentialFix::update()
 *          function is of no use here.
 * @see DataClass, likelihood_exponential, priormixlik_poisson,
 *      LogExponentialFix::update()
 * @author Lars Simon Zehnder
 * -------------------------------------------------------------
 **/
void LogExponentialInd::update(const unsigned int& K,
                               const arma::mat& y, arma::ivec &S, const arma::mat& expos,
                               const arma::vec& T, const ParExponentialInd& par,
                               const PriorExponentialInd& hyperPar)
{
   liklist   lik   = likelihood_exponential(y, par.lambda);
   DataClass dataC = classification(S, lik, par.weight);

   S      = dataC.newS;
   mixlik = dataC.mixLik;
   /* Compute likelihood of mixture prior */
   mixprior = priormixlik_poisson(par.lambda,
                                  hyperPar.aStart, hyperPar.bStart,
                                  hyperPar.HIER, hyperPar.g, hyperPar.G);
   if (K > 1)
   {
      /* Compute likelihood of Dirichlet prior */
      mixprior += priormixlik_dirichlet(par.weight,
                                        hyperPar.weightStart);
      cdpost  = mixlik + mixprior + dataC.postS;
      entropy = dataC.entropy;
   }
}
#endif
