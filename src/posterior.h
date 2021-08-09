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
#ifndef POSTERIOR_H
#define POSTERIOR_H

#include <RcppArmadillo.h>

/**
 * posterior for multinomial distribution
 *
 * as this is used for the weights of any mixture 
 * no specific struct is input argument 
 */

inline arma::rowvec 
posterior_multinomial(const unsigned int &K, const arma::ivec &S, 
			const arma::rowvec &weight) 
{
	arma::imat repS = arma::repmat(S, 1, K);
	arma::imat compM = arma::ones<arma::imat>(S.n_rows, K);
	arma::rowvec par_post(K);
 
	/* create sequence */
	for(unsigned int k = 0; k < K; ++k) {
		compM.col(k) = compM.col(k) * (k + 1);		
	}
	arma::umat ind = (repS == compM);
	arma::mat indDouble = arma::conv_to<arma::mat>::from(ind);
	par_post = arma::sum(indDouble); 
	par_post = par_post + weight;
	
	return par_post;
}

#endif
