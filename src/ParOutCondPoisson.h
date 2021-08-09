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
 * along with 'finmix'. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#ifndef __FINMIX_PAROUTCONDPOISSON_H_
#define __FINMIX_PAROUTCONDPOISSON_H_

#include <RcppArmadillo.h>
#include "ParCondPoissonFix.h"

class ParOutCondPoisson {
	public:
		arma::mat* lambda;
        Rcpp::NumericVector acc;

		ParOutCondPoisson () {}		
		ParOutCondPoisson (const Rcpp::List&); 
		void store(const unsigned int&, const ParCondPoissonFix&);
};

// =============================================================
// Constructor
// -------------------------------------------------------------

/**
 * -------------------------------------------------------------
 * ParOutPoisson
 * @brief   Constructs a ParOutiCondPoisson from an Rcpp::List
 *          object. Reuses memory allocated in R.
 * @par list    Rcpp::List object containing an R 'array'
 *              object, M x K, to store the sampled para-
 *              meters
 * @return  ParOutPoisson object
 * @detail  Reusage of memory allocated in R is done via the 
 *          Rcpp API and passing a pointer to the Armadillo
 *          matrix.
 * @see arma::mat::mat(), Rcpp::List
 * @author Lars Simon Zehnder
 * ------------------------------------------------------------
 **/
ParOutCondPoisson::ParOutCondPoisson(const Rcpp::List& list) :
    acc((SEXP) list["acc"])
{
	Rcpp::NumericMatrix tmpLambda((SEXP) list["lambda"]);
	const unsigned int M = tmpLambda.nrow();
	const unsigned int K = tmpLambda.ncol();
	lambda = new arma::mat(tmpLambda.begin(), M, K, false, true);      
} 

// =============================================================
// Store
// -------------------------------------------------------------

/** 
 * -------------------------------------------------------------
 * store
 * @brief   Stores the sampled parameters from step 'm'.
 * @par m   iteration step
 * @par par object of class ParCondPoissonFix holding the sampled
 *          parameters from step 'm'
 * @see ParPoissonFix
 * @author Lars Simon Zehnder
 * -------------------------------------------------------------
 **/
void ParOutCondPoisson::store (const unsigned int& m, const ParCondPoissonFix& par)
{
	(*lambda).row(m) = par.lambda;
    acc(0) = acc(0) + par.acc / (double) lambda->n_rows;
}

#endif
