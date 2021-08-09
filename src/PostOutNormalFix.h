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

#ifndef __FINMIX_POSTOUTNORMALFIX_H__
#define __FINMIX_POSTOUTNORMALFIX_H__

#include "PriorNormalFix.h"

class PostOutNormalFix {
    public:
            arma::mat* b;
            arma::mat* B;
            arma::mat* c;
            arma::mat* C;

            PostOutNormalFix () {}
            PostOutNormalFix (const Rcpp::List&);
            ~PostOutNormalFix () {}
            void store (const unsigned int&, 
                    const PriorNormalFix&);
};

PostOutNormalFix::PostOutNormalFix (const Rcpp::List& list)
{
    Rcpp::List tmpList((SEXP) list["par"]);
    Rcpp::List tmpMu((SEXP) tmpList["mu"]);
    Rcpp::List tmpSigma((SEXP) tmpList["sigma"]);
    Rcpp::NumericMatrix tmpb((SEXP) tmpMu["b"]);
    Rcpp::NumericMatrix tmpB((SEXP) tmpMu["B"]);
    Rcpp::NumericMatrix tmpc((SEXP) tmpSigma["c"]);
    Rcpp::NumericMatrix tmpC((SEXP) tmpSigma["C"]);
    const unsigned int M = tmpb.nrow();
    const unsigned int K = tmpb.ncol();
    b = new arma::mat(tmpb.begin(), M, K, false, true);
    B = new arma::mat(tmpB.begin(), M, K, false, true);
    c = new arma::mat(tmpc.begin(), M, K, false, true);
    C = new arma::mat(tmpC.begin(), M, K, false, true);
}

inline
void PostOutNormalFix::store (const unsigned int& m, 
        const PriorNormalFix& hyperPar)
{
    (*b).row(m) = hyperPar.bPost;
    (*B).row(m) = hyperPar.BPost;
    (*c).row(m) = hyperPar.cPost;
    (*C).row(m) = hyperPar.CPost;
}
#endif /* __FINMIX_POSTOUTNORMALFIX_H__ */



