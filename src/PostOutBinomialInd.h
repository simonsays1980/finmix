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

#ifndef __FINMIX_POSTOUTBINOMIALIND_H__
#define __FINMIX_POSTOUTBINOMIALIND_H__

#include "PostOutBinomialFix.h"
#include "PriorBinomialInd.h"

class PostOutBinomialInd : public PostOutBinomialFix {
    public:
        arma::mat* weight;

        PostOutBinomialInd (const Rcpp::List&);
        virtual ~PostOutBinomialInd () {}
        virtual void store (const unsigned int& m,
                const PriorBinomialInd&); 
            
};

PostOutBinomialInd::PostOutBinomialInd (const Rcpp::List& list) :
    PostOutBinomialFix(list) 
{
    Rcpp::NumericMatrix tmpWeight((SEXP) list["weight"]);
    const unsigned int M = tmpWeight.nrow();
    const unsigned int K = tmpWeight.ncol();
    weight = new arma::mat(tmpWeight.begin(), M, K, false, true);
}

void PostOutBinomialInd::store(const unsigned int& m,
        const PriorBinomialInd& hyperPar)
{
    PostOutBinomialFix::store(m, hyperPar);
    (*weight).row(m) = hyperPar.weightPost;
}
#endif /* __FINMIX_POSTOUTBINOMIALIND_H__ */



