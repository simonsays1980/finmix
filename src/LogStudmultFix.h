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

#ifndef __FINMIX_LOGSTUDMULTFIX_H__
#define __FINMIX_LOGSTUDMULTFIX_H__

#include "ParStudmultFix.h"
#include "likelihood.h"
#include "DataClass.h"
#include "prior_likelihood.h"

class LogStudmultFix {
    public:
        double mixlik;
        double mixprior;

        LogStudmultFix ();
        virtual ~LogStudmultFix () {}
        void update (const unsigned int&, const arma::mat&,
                const arma::ivec&, const arma::mat&,
                const arma::vec&, const ParStudmultFix&,
                const PriorStudmultFix&);
};

LogStudmultFix::LogStudmultFix () : mixlik(0.0),
    mixprior(0.0) {}

inline
void LogStudmultFix::update (const unsigned int& K, const arma::mat& y,
        const arma::ivec& S, const arma::mat& expos, 
        const arma::vec& T, const ParStudmultFix& par,
        const PriorStudmultFix& hyperPar)
{
    liklist lik     = likelihood_studmult(y, par.mu, par.sigma, par.df);
    DataClass dataC = classification_fix(K, S, lik);
    mixlik          = arma::sum(dataC.logLikCd);
    if (mixlik > 0.0) { 
        hyperPar.bPost.print("bPost");
        hyperPar.BPost.print("BPost");
        hyperPar.cPost.print("cPost");
        hyperPar.CPost.print("CPost");
    }
    /* Compute likelihood of mixture prior */
    mixprior        = priormixlik_studmult(hyperPar.INDEPENDENT, hyperPar.HIER,
            hyperPar.bStart, hyperPar.BInvStart, hyperPar.BStart, hyperPar.cStart, 
            hyperPar.CStart, hyperPar.logdetC, hyperPar.g, hyperPar.G, 
            par.mu, par.sigma, par.df, hyperPar.trans,
            hyperPar.a0, hyperPar.b0, hyperPar.d);

}
#endif /* __FINMIX_LOGSTUDENTFIX_H__ */



