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

#ifndef __FINMIX_PAROUTNORMAL_H__
#define __FINMIX_PAROUTNORMAL_H__

#include "ParNormalFix.h"

class ParOutNormal {
    public:
        arma::mat* mu;
        arma::mat* sigma;

        ParOutNormal () {}
        ParOutNormal (const Rcpp::List&);
        ~ParOutNormal () {}
        void store (const unsigned int&, const ParNormalFix&);
};

ParOutNormal::ParOutNormal (const Rcpp::List& list)
{
    Rcpp::NumericMatrix tmpMu((SEXP) list["mu"]);
    Rcpp::NumericMatrix tmpSigma((SEXP) list["sigma"]);
    const unsigned int M = tmpMu.nrow();
    const unsigned int K = tmpMu.ncol();
    mu      = new arma::mat(tmpMu.begin(), M, K, false, true);
    sigma   = new arma::mat(tmpSigma.begin(), M, K, false, true);    
}

void ParOutNormal::store (const unsigned int& m, const ParNormalFix& par)
{
    (*mu).row(m)    = par.mu;
    (*sigma).row(m) = par.sigma;
}
#endif /* __FINMIX_PAROUTNORMAL_H__ */



