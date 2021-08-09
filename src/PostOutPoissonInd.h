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
#ifndef __FINMIX_POSTOUTPOISSONIND_H__
#define __FINMIX_POSTOUTPOISSONIND_H__

#include "PostOutPoissonFix.h"
#include "PriorPoissonInd.h"

class PostOutPoissonInd : public PostOutPoissonFix {
	public:
		arma::mat* weight;
		
		PostOutPoissonInd (const Rcpp::List&);
		virtual ~PostOutPoissonInd () {}
		virtual void store (const unsigned int&,
				const PriorPoissonInd&); 
};

PostOutPoissonInd::PostOutPoissonInd (const Rcpp::List& list) :
	PostOutPoissonFix(list) 
{	
	Rcpp::NumericMatrix tmpWeight((SEXP) list["weight"]);
	const unsigned int M = tmpWeight.nrow();
	const unsigned int K = tmpWeight.ncol();
	weight = new arma::mat(tmpWeight.begin(), M, K, false, true);
}

void PostOutPoissonInd::store (const unsigned int& m,
	const PriorPoissonInd& hyperPar) 
{
	PostOutPoissonFix::store(m, hyperPar);
	(*weight).row(m) = hyperPar.weightPost;
}
#endif // __FINMIX_POSTOUTPOISSONIND_H_
