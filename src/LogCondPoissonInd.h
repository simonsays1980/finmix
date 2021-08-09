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
#ifndef __FINMIX_LOGCONDPOISSONIND_H_
#define __FINMIX_LOGCONDPOISSONIND_H_

#include <RcppArmadillo.h>
#include "LogCondPoissonFix.h"
#include "ParCondPoissonInd.h"

class LogCondPoissonInd : public LogCondPoissonFix {
	public:
		double cdpost;
		double entropy;
		double maxcdpost;

		LogCondPoissonInd ();
		virtual ~LogCondPoissonInd () {}
		void update (const unsigned int&, const arma::mat&,
			arma::ivec&, const arma::mat&, const arma::vec&,
            const ParCondPoissonInd&, const PriorCondPoissonInd&);
};

// =============================================================
// Constructor
// -------------------------------------------------------------
LogCondPoissonInd::LogCondPoissonInd () : LogCondPoissonFix(),
	cdpost(0.0), entropy(0.0), maxcdpost(0.0) {}

/** 
 * -------------------------------------------------------------
 * update
 * @brief   Updates the log-likelihoods of the Poisson model
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
 *          use of inheritance, i.e. the LogPoissonFix::update()
 *          function is of no use here.
 * @see DataClass, likelihood_poisson, priormixlik_poisson,
 *      LogPoissonFix::update()
 * @author Lars Simon Zehnder
 * -------------------------------------------------------------
 **/
void LogCondPoissonInd::update (const unsigned int& K, 
	const arma::mat& y, arma::ivec &S, const arma::mat& expos,
    const arma::vec& T, const ParCondPoissonInd& par, 
    const PriorCondPoissonInd& hyperPar)
{
	arma::mat lambdaM = arma::kron(expos, par.lambda);
    liklist lik = likelihood_poisson(y, lambdaM);
    DataClass dataC = classification(S, lik, par.weight);
    S = dataC.newS;
	mixlik = dataC.mixLik;
    /* Compute likelihood of mixture prior */
    mixprior = priormixlik_condpoisson(par.lambda, hyperPar);
    if(K > 1) {
		/* Compute likelihood of Dirichlet prior */
		mixprior += priormixlik_dirichlet(par.weight, 
			hyperPar.weightStart);
		cdpost = mixlik + mixprior + dataC.postS;
		entropy = dataC.entropy;
	}
	
}
#endif
