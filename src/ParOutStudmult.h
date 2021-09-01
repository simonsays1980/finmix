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

#ifndef __FINMIX_PAROUTSTUDMULT_H__
#define __FINMIX_PAROUTSTUDMULT_H__

#include "ParStudmultFix.h"
#include "mincol.h"

class ParOutStudmult {
    public:
        arma::cube*     mu;
        arma::cube*     sigma;
        arma::cube*     sigmainv;
        arma::mat*      df;
        arma::rowvec*   acc;
        unsigned int M;
        unsigned int r;
        unsigned int s;
        unsigned int K;
        bool STOREINV;

        ParOutStudmult () {}
        ParOutStudmult (const Rcpp::List&);
        ~ParOutStudmult () {}
        void store (const unsigned int&, const ParStudmultFix&);
};

ParOutStudmult::ParOutStudmult (const Rcpp::List& list) :
     M(0), r(0), s(0), K(0),STOREINV(false)
{
    STOREINV = Rcpp::as<bool>(list["storeinv"]);
    /* mu is an (M x r x K) array */ 
    Rcpp::NumericVector tmpMu((SEXP) list["mu"]);
    Rcpp::IntegerVector tmpMuDim = tmpMu.attr("dim");
    M = tmpMuDim[0];
    r = tmpMuDim[1];
    s = r * (r + 1) / 2;
    K = tmpMuDim[2];
    mu      = new arma::cube(tmpMu.begin(), M, r, K, false, true);    
    /* sigma is an (M x r(r + 1)/2 x K) array */
    Rcpp::NumericVector tmpSigma((SEXP) list["sigma"]);
    sigma   = new arma::cube(tmpSigma.begin(), M, s, K, false, true);
    if (STOREINV) {
        /* sigmainv is an (M x r(r + 1)/2 x K) array */
        Rcpp::NumericVector tmpSigmaInv((SEXP) list["sigmainv"]);        
        sigmainv = new arma::cube(tmpSigmaInv.begin(), M, s, K, false, true);
    }
    /* df is an M x K array */
    Rcpp::NumericMatrix tmpDf((SEXP) list["df"]);
    Rcpp::NumericVector tmpAcc((SEXP) list["acc"]);
    df = new arma::mat(tmpDf.begin(), M, K, false, true);
    acc = new arma::rowvec(tmpAcc.begin(), K, false, true);
}

void ParOutStudmult::store (const unsigned int& m, const ParStudmultFix& par)
{
    /* mu is a r x K matrix */
    mu->tube(m, 0, m, r - 1)    = par.mu;
    /* sigma is a cube and is transformed to an r * (r + 1) / 2  matrix */
    sigma->tube(m, 0, m, s - 1)   = cincolmat(par.sigma);
    if (STOREINV) {
        /* sigmainv is a cube and is transformed to an r * (r + 1) / 2 matrix */
        arma::mat tmp = cincolmat(par.sigmainv);
        sigmainv->tube(m, 0, m, s - 1) = cincolmat(par.sigmainv);
    }
    df->row(m) = par.df;
    *acc        = *acc + par.acc / (double) mu->n_rows;
}
#endif /* __FINMIX_PAROUTSTUDMULT_H__ */



