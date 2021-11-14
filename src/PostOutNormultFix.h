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

#ifndef __FINMIX_POSTOUTNORMULTFIX_H__
#define __FINMIX_POSTOUTNORMULTFIX_H__

#include "PriorNormultFix.h"
#include "mincol.h"

class PostOutNormultFix {
public:
arma::cube* b;
arma::cube* B;
arma::mat* c;
arma::cube* C;
unsigned int M;
unsigned int r;
unsigned int s;
unsigned int K;

PostOutNormultFix ()
{
}
PostOutNormultFix (const Rcpp::List&);
virtual ~PostOutNormultFix ()
{
}
void store(const unsigned int&,
           const PriorNormultFix&);
};

inline
PostOutNormultFix::PostOutNormultFix (const Rcpp::List& list) :
   M(0), r(0), s(0), K(0)
{
   Rcpp::List          tmpList((SEXP)list["par"]);
   Rcpp::List          tmpMu((SEXP)tmpList["mu"]);
   /* b is an (M x r x K) array */
   Rcpp::NumericVector tmpb((SEXP)tmpMu["b"]);
   Rcpp::IntegerVector tmpbDim = tmpb.attr("dim");

   M = tmpbDim[0];
   r = tmpbDim[1];
   s = r * (r + 1) / 2;
   K = tmpbDim[2];
   b = new arma::cube(tmpb.begin(), M, r, K, false, true);
   /* B is an (M x r(r + 1)/2 x K) array */
   Rcpp::NumericVector tmpB((SEXP)tmpMu["B"]);

   B = new arma::cube(tmpB.begin(), M, r * (r + 1) / 2, K, false, true);

   Rcpp::List          tmpSigma((SEXP)tmpList["sigma"]);
   /* c is an (M x K) array */
   Rcpp::NumericMatrix tmpc((SEXP)tmpSigma["c"]);

   c = new arma::mat(tmpc.begin(), M, K, false, true);
   /* C is an (M x r(r + 1)/2 x K) array */
   Rcpp::NumericVector tmpC((SEXP)tmpSigma["C"]);

   C = new arma::cube(tmpC.begin(), M, r * (r + 1) / 2, K, false, true);
}

inline
void PostOutNormultFix::store(const unsigned int& m,
                              const PriorNormultFix& hyperPar)
{
   b->tube(m, 0, m, r - 1) = hyperPar.bPost;
   B->tube(m, 0, m, s - 1) = cincolmat(hyperPar.BPost);
   c->row(m)               = hyperPar.cPost;
   C->tube(m, 0, m, s - 1) = cincolmat(hyperPar.CPost);
}
#endif /* __FINMIX_POSTOUTNORMULTFIX_H__ */



