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

#ifndef __FINMIX_LOGNORMULTFIX_H__
#define __FINMIX_LOGNORMULTFIX_H__

#include "ParNormultFix.h"
#include "likelihood.h"
#include "DataClass.h"
#include "prior_likelihood.h"

class LogNormultFix {
    public:
        double mixlik;
        double mixprior;

        LogNormultFix ();
        virtual ~LogNormultFix () {}
        void update (const unsigned int&, const arma::mat&,
                const arma::ivec&, const arma::mat&,
                const arma::vec&, const ParNormultFix&,
                const PriorNormultFix&);
};

inline
LogNormultFix::LogNormultFix () : mixlik(0.0),
    mixprior(0.0) {}

inline
void LogNormultFix::update (const unsigned int& K, const arma::mat& y,
        const arma::ivec& S, const arma::mat& expos,
        const arma::vec& T, const ParNormultFix& par,
        const PriorNormultFix& hyperPar) 
{
    liklist lik = likelihood_normult(y, par.mu, par.sigma);
    DataClass dataC = classification_fix(K, S, lik);
    mixlik = arma::sum(dataC.logLikCd);
    /* Compute likelihood of mixture prior */
    mixprior = priormixlik_normult(hyperPar.INDEPENDENT,
            hyperPar.HIER, hyperPar.bStart, hyperPar.BInvStart,
            hyperPar.BStart, hyperPar.cStart, hyperPar.CStart, 
            hyperPar.logdetC, hyperPar.g, hyperPar.G, par.mu, 
            par.sigma);
}
#endif /* __FINMIX_LOGNORMULTFIX_H__ */



