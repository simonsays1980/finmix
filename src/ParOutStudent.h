/******************************************************************************
*
* TODO: Project Title
*
* Copyright (C) 2003-2009 ascolab GmbH. All Rights Reserved.
* Web: http://www.ascolab.com
*
* Author: Gerhard Gappmeier <gerhard.gappmeier@ascolab.com>
*
* This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
* WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
*
******************************************************************************/

#ifndef __FINMIX_PAROUTSTUDENT_H__
#define __FINMIX_PAROUTSTUDENT_H__

#include "ParStudentFix.h"

class ParOutStudent {
public:
arma::mat* mu;
arma::mat* sigma;
arma::mat* df;
arma::rowvec* acc;

ParOutStudent ()
{
}
ParOutStudent (const Rcpp::List& list);
~ParOutStudent ()
{
}
void store(const unsigned int&, const ParStudentFix&);
};

ParOutStudent::ParOutStudent (const Rcpp::List& list)
{
   Rcpp::NumericMatrix tmpMu((SEXP)list["mu"]);
   Rcpp::NumericMatrix tmpSigma((SEXP)list["sigma"]);
   Rcpp::NumericMatrix tmpDf((SEXP)list["df"]);
   Rcpp::NumericVector tmpAcc((SEXP)list["acc"]);
   const unsigned int  M = tmpMu.nrow();
   const unsigned int  K = tmpMu.ncol();

   mu    = new arma::mat(tmpMu.begin(), M, K, false, true);
   sigma = new arma::mat(tmpSigma.begin(), M, K, false, true);
   df    = new arma::mat(tmpDf.begin(), M, K, false, true);
   acc   = new arma::rowvec(tmpAcc.begin(), K, false, true);
}

inline
void ParOutStudent::store(const unsigned int& m, const ParStudentFix& par)
{
   (*mu).row(m)    = par.mu;
   (*sigma).row(m) = par.sigma;
   (*df).row(m)    = par.df;
   *acc            = *acc + par.acc / (double)mu->n_rows;
}
#endif /* __FINMIX_PAROUTSTUDENT_H__ */



