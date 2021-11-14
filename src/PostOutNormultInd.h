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

#ifndef __FINMIX_POSTOUTNORMULTIND_H__
#define __FINMIX_POSTOUTNORMULTIND_H__

#include "PostOutNormultFix.h"
#include "PriorNormultInd.h"

class PostOutNormultInd : public PostOutNormultFix {
public:
arma::mat* weight;

PostOutNormultInd (const Rcpp::List&);
virtual ~PostOutNormultInd ()
{
}
using PostOutNormultFix::store;
void store(const unsigned int&,
           const PriorNormultInd&);
};

inline
PostOutNormultInd::PostOutNormultInd (const Rcpp::List& list) :
   PostOutNormultFix(list)
{
   Rcpp::NumericMatrix tmpWeight((SEXP)list["weight"]);

   weight = new arma::mat(tmpWeight.begin(), M, K, false, true);
}

inline
void PostOutNormultInd::store(const unsigned int& m,
                              const PriorNormultInd& hyperPar)
{
   PostOutNormultFix::store(m, hyperPar);
   (*weight).row(m) = hyperPar.weightPost;
}
#endif /* __FINMIX_POSTOUTNORMULTIND_H__ */



