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
#ifndef __FINMIX_HIEROUTEXPONENTIAL_H__
#define __FINMIX_HIEROUTEXPONENTIAL_H__

#include <RcppArmadillo.h>

class HierOutExponential {
public:
arma::vec* b;

HierOutExponential ()
{
}
HierOutExponential (const Rcpp::List&);
template <typename PriorParType>
void store(const unsigned int& m,
           const PriorParType&);
};

inline
HierOutExponential::HierOutExponential (const Rcpp::List& list)
{
   Rcpp::NumericVector tmpB((SEXP)list["b"]);
   const unsigned int  M = tmpB.size();

   b = new arma::vec(tmpB.begin(), M, false, true);
}

template <typename PriorParType>
inline
void HierOutExponential::store(const unsigned int& m,
                               const PriorParType& hyperPar)
{
   (*b)(m) = hyperPar.bStart(0);
}
#endif // __FINMIX_HIEROUTEXPONENTIAL_H__
