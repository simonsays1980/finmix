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

#ifndef __FINMIX_LOGNORMALFIX_H__
#define __FINMIX_LOGNORMALFIX_H__

#include "ParNormalFix.h"
#include "likelihood.h"
#include "DataClass.h"
#include "prior_likelihood.h"

class LogNormalFix {
    public:
        double mixlik;
        double mixprior;

        LogNormalFix ();
        virtual ~LogNormalFix () {}
        void update (const unsigned int&, const arma::mat&,
                const arma::ivec&, const arma::mat&,
                const arma::vec&, const ParNormalFix&,
                const PriorNormalFix&);
};

LogNormalFix::LogNormalFix () : mixlik(0.0), 
    mixprior(0.0) {}

void LogNormalFix::update (const unsigned int& K, const arma::mat& y,
        const arma::ivec& S, const arma::mat& expos, 
        const arma::vec& T, const ParNormalFix& par,
        const PriorNormalFix& hyperPar)
{
    liklist lik = likelihood_normal(y, par.mu, par.sigma);
    DataClass dataC = classification_fix(K, S, lik);
    mixlik = arma::sum(dataC.logLikCd);
    /* Compute likelihood of mixture prior */
    mixprior = priormixlik_normal(hyperPar.INDEPENDENT,
            hyperPar.HIER,
            hyperPar.bStart, hyperPar.BStart, 
            hyperPar.cStart, hyperPar.CStart,
            par.mu, par.sigma, hyperPar.g, hyperPar.G);
}
#endif /* __FINMIX_LOGNORMALFIX_H__ */



