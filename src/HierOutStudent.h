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

#ifndef __FINMIX_HIEROUTSTUDENT_H__
#define __FINMIX_HIEROUTSTUDENT_H__

class HierOutStudent {
    public:
        arma::vec* C;

        HierOutStudent () {}
        HierOutStudent (const Rcpp::List&);
        ~HierOutStudent () {}
        template <typename PriorParType>
        void store (const unsigned int&, 
                const PriorParType&);
};

HierOutStudent::HierOutStudent (const Rcpp::List& list)
{
    Rcpp::NumericVector tmpC((SEXP) list["C"]);
    const unsigned int M = tmpC.size();
    C = new arma::vec(tmpC.begin(), M, false, true);
}

template <typename PriorParType>
inline
void HierOutStudent::store (const unsigned int& m,
        const PriorParType& hyperPar)
{
    (*C)(m) = hyperPar.CStart(0);
}

#endif /* __FINMIX_HIEROUTSTUDENT_H__ */



