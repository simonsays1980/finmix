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
#include "PriorExponentialFix.h"
#include "ParExponentialFix.h"

PriorExponentialFix::PriorExponentialFix () : HIER(false)
{
}


PriorExponentialFix::PriorExponentialFix (const FinmixPrior& prior) :
   aStart(Rcpp::as<arma::rowvec>((SEXP)prior.par["a"])),
   bStart(Rcpp::as<arma::rowvec>((SEXP)prior.par["b"])),
   aPost(Rcpp::as<arma::rowvec>((SEXP)prior.par["a"])),
   bPost(Rcpp::as<arma::rowvec>((SEXP)prior.par["b"])),
   HIER(false)
{
}

void PriorExponentialFix::update(const unsigned int& K, const arma::mat& y,
                                 arma::ivec& S, const arma::vec& T, const ParExponentialFix& par)
{
   if (K == 1)
   {
      aPost(0) = aStart(0) + arma::accu(y);
      bPost(0) = bStart(0) + y.n_rows;
   }
   else
   {
      arma::mat  repY  = arma::repmat(y, 1, K);
      arma::imat repS  = arma::repmat(S, 1, K);
      arma::imat compM = arma::ones<arma::imat>(S.n_elem, K);
      for (unsigned int k = 0; k < K; ++k)
      {
         compM.col(k) = compM.col(k) * (k + 1);
      }
      arma::umat   ind       = (repS == compM);
      arma::mat    indDouble = arma::conv_to<arma::mat>::from(ind);
      repY %= indDouble;
      arma::rowvec sprod = sum(repY, 0);
      arma::rowvec sind  = sum(indDouble, 0);
      aPost = aStart + sind;
      bPost = bStart + sprod;
   }
}

void PriorExponentialFix::updateHier(const ParExponentialFix& par)
{
}

