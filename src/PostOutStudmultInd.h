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

#ifndef __FINMIX_POSTOUTSTUDMULTIND_H__
#define __FINMIX_POSTOUTSTUDMULTIND_H__

#include "PostOutStudmultFix.h"
#include "PriorStudmultInd.h"

class PostOutStudmultInd : public PostOutStudmultFix {
public:
arma::mat* weight;

PostOutStudmultInd (const Rcpp::List&);
virtual ~PostOutStudmultInd ()
{
}
virtual void store(const unsigned int&,
                   const PriorStudmultInd&);
};

inline
PostOutStudmultInd::PostOutStudmultInd (const Rcpp::List& list) :
   PostOutStudmultFix(list)
{
   Rcpp::NumericMatrix tmpWeight((SEXP)list["weight"]);

   weight = new arma::mat(tmpWeight.begin(), M, K, false, true);
}

inline
void PostOutStudmultInd::store(const unsigned int& m,
                               const PriorStudmultInd& hyperPar)
{
   PostOutStudmultFix::store(m, hyperPar);
   (*weight).row(m) = hyperPar.weightPost;
}
#endif /* __FINMIX_POSTOUTSTUDMULTIND_H__ */



