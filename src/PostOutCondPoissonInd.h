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
#ifndef __FINMIX_POSTOUTCONDPOISSONIND_H_
#define __FINMIX_POSTOUTCONDPOISSONIND_H_

#include <RcppArmadillo.h>
#include "PostOutCondPoissonFix.h"
#include "PriorCondPoissonInd.h"

class PostOutCondPoissonInd : public PostOutCondPoissonFix {
	public:
		arma::mat* weight;
		
		PostOutCondPoissonInd (const Rcpp::List&);
		virtual ~PostOutCondPoissonInd () {}
		virtual void store (const unsigned int&,
				const PriorCondPoissonInd&); 
};

PostOutCondPoissonInd::PostOutCondPoissonInd (const Rcpp::List& list) :
	PostOutCondPoissonFix(list) 
{	
	Rcpp::NumericMatrix tmpWeight((SEXP) list["weight"]);
	const unsigned int M = tmpWeight.nrow();
	const unsigned int K = tmpWeight.ncol();
	weight = new arma::mat(tmpWeight.begin(), M, K, false, true);
}

void PostOutCondPoissonInd::store (const unsigned int& m,
	const PriorCondPoissonInd& hyperPar) 
{
	PostOutCondPoissonFix::store(m, hyperPar);
	(*weight).row(m) = hyperPar.weightPost;
}
#endif // __FINMIX_POSTOUTCONDPOISSOND_H_
