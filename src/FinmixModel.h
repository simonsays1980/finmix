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
#ifndef FINMIXMODEL_H
#define FINMIXMODEL_H

#include <RcppArmadillo.h>
#include <string>

class FinmixModel {
public:
Rcpp::List par;
arma::rowvec weight;

bool indicFix;
unsigned int K;
unsigned int r;

/* ctor */
FinmixModel (Rcpp::S4& classS4) :
   par(classS4.slot("par")),
   weight(Rcpp::as<arma::mat>((SEXP)classS4.slot("weight"))),
   indicFix(Rcpp::as<bool>((SEXP)classS4.slot("indicfix"))),
   K(Rcpp::as<unsigned int>((SEXP)classS4.slot("K"))),
   r(Rcpp::as<unsigned int>((SEXP)classS4.slot("r")))
{
}

/* dtor */
~FinmixModel()
{
}
};
#endif
