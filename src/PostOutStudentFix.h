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

#ifndef __FINMIX_POSTOUTSTUDENTFIX_H__
#define __FINMIX_POSTOUTSTUDENTFIX_H__

#include "PriorStudentFix.h"

class PostOutStudentFix {
public:
arma::mat* b;
arma::mat* B;
arma::mat* c;
arma::mat* C;

PostOutStudentFix ()
{
}
PostOutStudentFix (const Rcpp::List&);
~PostOutStudentFix ()
{
}
void store(const unsigned int&,
           const PriorStudentFix&);
};

PostOutStudentFix::PostOutStudentFix (const Rcpp::List& list)
{
   Rprintf("PostOut, Line 36\n");
   Rcpp::List          tmpList((SEXP)list["par"]);
   Rcpp::List          tmpMu((SEXP)tmpList["mu"]);
   Rcpp::List          tmpSigma((SEXP)tmpList["sigma"]);
   Rcpp::NumericMatrix tmpb((SEXP)tmpMu["b"]);
   Rcpp::NumericMatrix tmpB((SEXP)tmpMu["B"]);
   Rcpp::NumericMatrix tmpc((SEXP)tmpSigma["c"]);
   Rcpp::NumericMatrix tmpC((SEXP)tmpSigma["C"]);

   Rprintf("PostOut, Line 44\n");
   const unsigned int M = tmpb.nrow();
   const unsigned int K = tmpb.ncol();

   b = new arma::mat(tmpb.begin(), M, K, false, true);
   B = new arma::mat(tmpB.begin(), M, K, false, true);
   c = new arma::mat(tmpc.begin(), M, K, false, true);
   C = new arma::mat(tmpC.begin(), M, K, false, true);
}

inline
void PostOutStudentFix::store(const unsigned int& m,
                              const PriorStudentFix& hyperPar)
{
   (*b).row(m) = hyperPar.bPost;
   (*B).row(m) = hyperPar.BPost;
   (*c).row(m) = hyperPar.cPost;
   (*C).row(m) = hyperPar.CPost;
}

#endif /* __FINMIX_POSTOUTSTUDENTFIX_H__ */



