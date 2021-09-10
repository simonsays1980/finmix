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

#ifndef __FINMIX_HIEROUTNORMULT_H__
#define __FINMIX_HIEROUTNORMULT_H__

#include "mincol.h"

class HierOutNormult {
public:
arma::mat* C;
unsigned int s;

HierOutNormult ()
{
}
HierOutNormult (const Rcpp::List&);
~HierOutNormult ()
{
}
template <typename PriorParType>
void store(const unsigned int&,
           const PriorParType&);
};

inline
HierOutNormult::HierOutNormult (const Rcpp::List& list) :
   s(0)
{
   Rcpp::NumericMatrix tmpC((SEXP)list["C"]);

   s = tmpC.ncol();
   C = new arma::mat(tmpC.begin(), tmpC.nrow(), s, false, true);
}

template <typename PriorParType>
inline
void HierOutNormult::store(const unsigned int& m,
                           const PriorParType& hyperPar)
{
   C->row(m) = minrow(hyperPar.CStart.slice(0));
}
#endif /* __FINMIX_HIEROUTNORMULT_H__ */



