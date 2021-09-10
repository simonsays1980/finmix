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
#ifndef __FINMIX_PAROUTBINOMIAL_H__
#define __FINMIX_PAROUTBINOMIAL_H__

#include "ParBinomialFix.h"

class ParOutBinomial {
public:
arma::mat* p;

ParOutBinomial ()
{
}
ParOutBinomial (const Rcpp::List&);
void store(const unsigned int&, const ParBinomialFix&);
};

// =============================================================
// Constructor
// -------------------------------------------------------------

/**
 * -------------------------------------------------------------
 * ParOutBinomial
 * @brief   Constructs a ParOutBinomial from an Rcpp::List
 *          object. Reuses memory allocated in R.
 * @par list    Rcpp::List object containing an R 'array'
 *              object, M x K, to store the sampled para-
 *              meters
 * @return  ParOutBinomial object
 * @detail  reusage of memory allocated in R is done via the
 *          Rcpp API and passing apointer to the Armadillo
 *          matrix
 * @see arma::mat::mat(), Rcpp::List
 * @author Lars Simon Zehnder
 * ------------------------------------------------------------
 **/
ParOutBinomial::ParOutBinomial (const Rcpp::List& list)
{
   Rcpp::NumericMatrix tmpP((SEXP)list["p"]);
   const unsigned int  M = tmpP.nrow();
   const unsigned int  K = tmpP.ncol();

   p = new arma::mat(tmpP.begin(), M, K, false, true);
}

// =============================================================
// Store
// -------------------------------------------------------------

/**
 * -------------------------------------------------------------
 * store
 * @brief   Stores the sampled parameters from step 'm'.
 * @par m   iteration step
 * @par par object of class ParBinomialFix holding the sampled
 *          parameters from step 'm'
 * @see ParBinomialFix
 * @author Lars Simon Zehnder
 * -------------------------------------------------------------
 **/
void ParOutBinomial::store(const unsigned int& m,
                           const ParBinomialFix& par)
{
   (*p).row(m) = par.p;
}
#endif /* __FINMIX_PAROUTBINOMIAL_H__ */



