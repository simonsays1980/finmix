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
#ifndef __FINMIX_POSTOUTBINOMIALFIX_H__
#define __FINMIX_POSTOUTBINOMIALFIX_H__

#include "PriorBinomialFix.h"

class PostOutBinomialFix {
public:
arma::mat *a;
arma::mat *b;

PostOutBinomialFix ()
{
}
PostOutBinomialFix (const Rcpp::List&);
~PostOutBinomialFix ()
{
}
void store(const unsigned int&,
           const PriorBinomialFix&);
};

// =============================================================
// Constructor
// -------------------------------------------------------------

/**
 * -------------------------------------------------------------
 * PostOutBinomialFix
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
PostOutBinomialFix::PostOutBinomialFix (const Rcpp::List& list)
{
   Rcpp::List          tmpList((SEXP)list["par"]);
   Rcpp::NumericMatrix tmpA((SEXP)tmpList["a"]);
   Rcpp::NumericMatrix tmpB((SEXP)tmpList["b"]);
   const unsigned int  M = tmpA.nrow();
   const unsigned int  K = tmpA.ncol();

   a = new arma::mat(tmpA.begin(), M, K, false, true);
   b = new arma::mat(tmpB.begin(), M, K, false, true);
}

// =============================================================
// Store
// -------------------------------------------------------------

/**
 * -------------------------------------------------------------
 * store
 * @brief   Stores the posterior hyper parameters from step 'm'.
 * @par m   iteration step
 * @par par object of class PriorBinomialFix holding the posterior
 *          hyper parameters from step 'm'
 * @see PriorBinomialFix
 * @author Lars Simon Zehnder
 * -------------------------------------------------------------
 **/
void PostOutBinomialFix::store(const unsigned int& m,
                               const PriorBinomialFix& hyperPar)
{
   (*a).row(m) = hyperPar.aPost;
   (*b).row(m) = hyperPar.bPost;
}
#endif /* __FINMIX_POSTOUTBINOMIALFIX_H__ */



