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

#ifndef __FINMIX_LOGSTUDENTFIX_H__
#define __FINMIX_LOGSTUDENTFIX_H__

#include "ParStudentFix.h"
#include "likelihood.h"
#include "DataClass.h"
#include "prior_likelihood.h"

class LogStudentFix {
    public:
        double mixlik;
        double mixprior;

        LogStudentFix ();
        virtual ~LogStudentFix () {}
        void update (const unsigned int&, const arma::mat&,
                const arma::ivec&, const arma::mat&,
                const arma::vec&, const ParStudentFix&,
                const PriorStudentFix&);
};

LogStudentFix::LogStudentFix () : mixlik(0.0),
    mixprior(0.0) {}

inline
void LogStudentFix::update (const unsigned int& K, const arma::mat& y,
        const arma::ivec& S, const arma::mat& expos, 
        const arma::vec& T, const ParStudentFix& par,
        const PriorStudentFix& hyperPar)
{
    liklist lik     = likelihood_student(y, par.mu, par.sigma, par.df);
    DataClass dataC = classification_fix(K, S, lik);
    mixlik          = arma::sum(dataC.logLikCd);
    /* Compute likelihood of mixture prior */
    mixprior        = priormixlik_student(hyperPar.INDEPENDENT, hyperPar.HIER,
            hyperPar.bStart, hyperPar.BStart, hyperPar.cStart, hyperPar.CStart,
            par.mu, par.sigma, hyperPar.g, hyperPar.G, par.df, hyperPar.trans,
            hyperPar.a0, hyperPar.b0, hyperPar.d);

}
#endif /* __FINMIX_LOGSTUDENTFIX_H__ */



