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

#ifndef __FINMIX_POSTOUTNORMALIND_H__
#define __FINMIX_POSTOUTNORMALIND_H__

#include "PostOutNormalFix.h"
#include "PriorNormalInd.h"

class PostOutNormalInd: public PostOutNormalFix {
    public: 
        arma::mat* weight;

        PostOutNormalInd (const Rcpp::List&);
        virtual ~PostOutNormalInd () {}
        virtual void store (const unsigned int&, 
                const PriorNormalInd&);
};

inline
PostOutNormalInd::PostOutNormalInd (const Rcpp::List& list) :
    PostOutNormalFix(list) 
{
    Rcpp::NumericMatrix tmpWeight((SEXP) list["weight"]);
    const unsigned int M = tmpWeight.nrow();
    const unsigned int K = tmpWeight.ncol();
    weight = new arma::mat(tmpWeight.begin(), M, K, false, true);
}

inline
void PostOutNormalInd::store (const unsigned int& m, 
        const PriorNormalInd& hyperPar)
{
    PostOutNormalFix::store(m, hyperPar);
    (*weight).row(m) = hyperPar.weightPost;
}
#endif /* __FINMIX_POSTOUTNORMALIND_H__ */



