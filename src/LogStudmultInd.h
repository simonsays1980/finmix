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

#ifndef __FINMIX_LOGSTUDMULTIND_H__
#define __FINMIX_LOGSTUDMULTIND_H__

#include "LogStudmultFix.h"
#include "ParStudmultFix.h"
#include "PriorStudmultFix.h"

class LogStudmultInd : public LogStudmultFix {
    public:
        double cdpost;
        double entropy;
        double maxcdpost;

        LogStudmultInd ();
        virtual ~LogStudmultInd () {}
        void update (const unsigned int&, const arma::mat&, 
                arma::ivec&, const arma::mat&, const arma::vec&,
                const ParStudmultInd&, const PriorStudmultInd&);
};

inline
LogStudmultInd::LogStudmultInd () : LogStudmultFix(),
    cdpost(0.0), entropy(0.0), maxcdpost(0.0) {}

inline
void LogStudmultInd::update (const unsigned int& K,
        const arma::mat& y, arma::ivec& S, const arma::mat& expos,
        const arma::vec& T, const ParStudmultInd& par,
        const PriorStudmultInd& hyperPar)
{
    liklist lik = likelihood_studmult(y, par.mu, par.sigmainv, par.df);
    DataClass dataC = classification(S, lik, par.weight);
    S = dataC.newS;
    mixlik = dataC.mixLik;

    /* Compute likelihood of mixture prior */
    mixprior        = priormixlik_studmult(hyperPar.INDEPENDENT, hyperPar.HIER,
            hyperPar.bStart, hyperPar.BInvStart, hyperPar.BStart, hyperPar.cStart, 
            hyperPar.CStart, hyperPar.logdetC, hyperPar.g, hyperPar.G, 
            par.mu, par.sigma, par.df, hyperPar.trans,
            hyperPar.a0, hyperPar.b0, hyperPar.d);
    if (K > 1) {
        /* Compute likelihood of DIrichlet prior */
        mixprior += priormixlik_dirichlet(par.weight, 
                hyperPar.weightStart);
        cdpost = mixlik + mixprior + dataC.postS;
        entropy = dataC.entropy;
    }
}
#endif /* __FINMIX_LOGSTUDMULTIND_H__ */



