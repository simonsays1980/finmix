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
#ifndef __FINMIX_FINMIXPRIOR_H__
#define __FINMIX_FINMIXPRIOR_H__

#include <RcppArmadillo.h>
#include <string>

class FinmixPrior {
public:
Rcpp::List par;
arma::rowvec weight;

std::string type;
bool hier;

/* ctor */
FinmixPrior (Rcpp::S4& classS4) :
   par(Rcpp::as<Rcpp::List>((SEXP)classS4.slot("par"))),
   weight(Rcpp::as<arma::rowvec>((SEXP)classS4.slot("weight"))),
   type(Rcpp::as<std::string>((SEXP)classS4.slot("type"))),
   hier(Rcpp::as<bool>((SEXP)classS4.slot("hier")))
{
}
};
#endif // __FINMIX_FINMIXPRIOR_H__
