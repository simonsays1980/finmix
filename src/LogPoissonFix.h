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
#ifndef LOGPOISSONFIX_H
#define LOGPOISSONFIX_H

#include <RcppArmadillo.h>
#include "ParPoissonFix.h"
#include "likelihood.h"
#include "DataClass.h"
#include "prior_likelihood.h"

class LogPoissonFix {
	public:
		double mixlik;
		double mixprior;
		
		LogPoissonFix ();
		virtual ~LogPoissonFix () {}
		void update (const unsigned int&, const arma::mat&, 
				const arma::ivec&, const arma::mat& expos,
                const arma::vec&, const ParPoissonFix&, 
                const PriorPoissonFix&);
};

LogPoissonFix::LogPoissonFix () : mixlik(0.0), 
				mixprior(0.0) {}

void LogPoissonFix::update (const unsigned int& K, const arma::mat& y, 
			const arma::ivec& S, const arma::mat& expos, 
            const arma::vec& T, const ParPoissonFix& par, 
			const PriorPoissonFix& hyperPar) 
{
	arma::mat lambdaM = arma::kron(expos, par.lambda);
	liklist lik = likelihood_poisson(y, lambdaM);
	DataClass dataC = classification_fix(K, S, lik);
	mixlik = arma::sum(dataC.logLikCd);
	/* Compute likelihood of mixture prior */
	mixprior = priormixlik_poisson(par.lambda, hyperPar.aStart, hyperPar.bStart,
			hyperPar.HIER, hyperPar.g, hyperPar.G);
}
#endif
