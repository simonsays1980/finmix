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
#include "PriorPoissonFix.h"
#include "ParPoissonFix.h"

PriorPoissonFix::PriorPoissonFix () : HIER(false) {}


PriorPoissonFix::PriorPoissonFix (const FinmixPrior& prior) :
	aStart(Rcpp::as<arma::rowvec>((SEXP) prior.par["a"])),
	bStart(Rcpp::as<arma::rowvec>((SEXP) prior.par["b"])),
	aPost(Rcpp::as<arma::rowvec>((SEXP) prior.par["a"])),
	bPost(Rcpp::as<arma::rowvec>((SEXP) prior.par["b"])),
	HIER(prior.hier) 
{
	if (HIER){
		g = Rcpp::as<double>((SEXP) prior.par["g"]);
		G = Rcpp::as<double>((SEXP) prior.par["G"]); 
	}
}

void PriorPoissonFix::update (const unsigned int& K, const arma::mat& y,
			arma::ivec& S, const arma::vec& T, const ParPoissonFix& par)  
{
	if (K == 1) {
		aPost(0) = aStart(0) + arma::accu(y);
		bPost(0) = bStart(0) + y.n_rows; 
	}
	else {
		arma::mat repY = arma::repmat(y, 1, K);
		arma::imat repS = arma::repmat(S, 1, K);
		arma::imat compM = arma::ones<arma::imat>(S.n_elem, K);
		for(unsigned int k = 0; k < K; ++k) {
			compM.col(k) = compM.col(k) * (k + 1);
		}
		arma::umat ind = (repS == compM);
		arma::mat indDouble = arma::conv_to<arma::mat>::from(ind);
		repY               %= indDouble;
		arma::rowvec sprod  = sum(repY, 0);
		arma::rowvec sind   = sum(indDouble, 0);
		aPost               = aStart + sprod;
		bPost               = bStart + sind;
	}
}

void PriorPoissonFix::updateHier(const ParPoissonFix& par)
{
	if(HIER){ // Hierarchical prior for 'b'
		GetRNGstate(); // Get RNG state from R
		// Sample from G(g_0 + Ka_0, G_0 + sum lambda_k)
		double gN = g + arma::sum(aStart);
		double GN = G + arma::sum(par.lambda);
		double b = R::rgamma(gN, 1/GN);
		PutRNGstate();
		bStart.fill(b);
	}
}

