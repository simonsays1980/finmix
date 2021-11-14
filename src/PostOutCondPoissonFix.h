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
#ifndef __FINMIX_POSTOUTPOISSONFIX_H_
#define __FINMIX_POSTOUTPOISSONFIX_H_

#include <RcppArmadillo.h>
#include "PriorCondPoissonFix.h"

class PostOutCondPoissonFix {
public:
arma::mat *Q;
arma::mat *N;

PostOutCondPoissonFix ()
{
}
PostOutCondPoissonFix (const Rcpp::List&);
~PostOutCondPoissonFix ()
{
}
void store(const unsigned int&,
           const PriorCondPoissonFix&);
};

PostOutCondPoissonFix::PostOutCondPoissonFix (const Rcpp::List& list)
{
   Rcpp::List          tmpList((SEXP)list["par"]);
   Rcpp::NumericMatrix tmpQ((SEXP)tmpList["Q"]);
   Rcpp::NumericMatrix tmpN((SEXP)tmpList["N"]);
   const unsigned int  M = tmpQ.nrow();
   const unsigned int  K = tmpQ.ncol();

   Q = new arma::mat(tmpQ.begin(), M, K, false, true);
   N = new arma::mat(tmpN.begin(), M, K, false, true);
}

void PostOutCondPoissonFix::store(const unsigned int& m,
                                  const PriorCondPoissonFix& hyperPar)
{
   (*Q).row(m) = hyperPar.Q;
   (*N).row(m) = hyperPar.N;
}
#endif
