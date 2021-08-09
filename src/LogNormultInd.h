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

#ifndef __FINMIX_LOGNORMULTIND_H__
#define __FINMIX_LOGNORMULTIND_H__

#include "LogNormultFix.h"
#include "ParNormultFix.h"
#include "PriorNormultFix.h"

class LogNormultInd : public LogNormultFix {
    public:
        double cdpost;
        double entropy;
        double maxcdpost;

        LogNormultInd ();
        virtual ~LogNormultInd () {}
        void update (const unsigned int&, const arma::mat&, 
                arma::ivec&, const arma::mat&, const arma::vec&,
                const ParNormultInd&, const PriorNormultInd&);
};

inline
LogNormultInd::LogNormultInd () : LogNormultFix(),
    cdpost(0.0), entropy(0.0), maxcdpost(0.0) {}

inline
void LogNormultInd::update (const unsigned int& K,
        const arma::mat& y, arma::ivec& S, const arma::mat& expos,
        const arma::vec& T, const ParNormultInd& par,
        const PriorNormultInd& hyperPar)
{
    liklist lik = likelihood_normult(y, par.mu, par.sigma);
    DataClass dataC = classification(S, lik, par.weight);
    S = dataC.newS;
    mixlik = dataC.mixLik;
    /* Compute likelihood of mixture prior */
    mixprior = priormixlik_normult(hyperPar.INDEPENDENT,
            hyperPar.HIER, hyperPar.bStart, hyperPar.BInvStart,
            hyperPar.BStart, hyperPar.cStart, hyperPar.CStart, 
            hyperPar.logdetC, hyperPar.g, hyperPar.G, par.mu,
            par.sigma);
    if (K > 1) {
        /* Compute likelihood of DIrichlet prior */
        mixprior += priormixlik_dirichlet(par.weight, 
                hyperPar.weightStart);
        cdpost = mixlik + mixprior + dataC.postS;
        entropy = dataC.entropy;
    }
}
#endif /* __FINMIX_LOGNORMULTIND_H__ */



