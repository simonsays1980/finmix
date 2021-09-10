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

#ifndef __FINMIX_POSTOUTSTUDENTIND_H__
#define __FINMIX_POSTOUTSTUDENTIND_H__

#include "PostOutStudentFix.h"
#include "PriorStudentInd.h"

class PostOutStudentInd : virtual public PostOutStudentFix {
public:
arma::mat* weight;

PostOutStudentInd (const Rcpp::List&);
virtual ~PostOutStudentInd ()
{
}
virtual void store(const unsigned int&,
                   const PriorStudentInd&);
};

PostOutStudentInd::PostOutStudentInd (const Rcpp::List& list) :
   PostOutStudentFix(list)
{
   Rcpp::NumericMatrix tmpWeight((SEXP)list["weight"]);
   const unsigned int  M = tmpWeight.nrow();
   const unsigned int  K = tmpWeight.ncol();

   weight = new arma::mat(tmpWeight.begin(), M, K, false, true);
}

void PostOutStudentInd::store(const unsigned int& m,
                              const PriorStudentInd& hyperPar)
{
   PostOutStudentFix::store(m, hyperPar);
   (*weight).row(m) = hyperPar.weightPost;
}

#endif /* __FINMIX_POSTOUTSTUDENTIND_H__ */



