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

#ifndef __FINMIX_LOGSTUDENTIND_H__
#define __FINMIX_LOGSTUDENTIND_H__

#include "LogStudentFix.h"
#include "ParNormalInd.h"
#include "PriorStudentInd.h"

class LogStudentInd : public LogStudentFix {
public:
double cdpost;
double entropy;
double maxcdpost;

LogStudentInd ();
virtual ~LogStudentInd ()
{
}
void update(const unsigned int&, const arma::mat&,
            arma::ivec&, const arma::mat&, const arma::vec&,
            const ParStudentInd&, const PriorStudentInd&);
};

LogStudentInd::LogStudentInd () : LogStudentFix(),
   cdpost(0.0), entropy(0.0), maxcdpost(0.0)
{
}

inline
void LogStudentInd::update(const unsigned int& K,
                           const arma::mat& y, arma::ivec& S, const arma::mat& expos,
                           const arma::vec& T, const ParStudentInd& par,
                           const PriorStudentInd& hyperPar)
{
   liklist   lik   = likelihood_student(y, par.mu, par.sigma, par.df);
   DataClass dataC = classification(S, lik, par.weight);

   S      = dataC.newS;
   mixlik = dataC.mixLik;
   /* Compute likelihood of mixture prior */
   mixprior = priormixlik_student(hyperPar.INDEPENDENT,
                                  hyperPar.HIER, hyperPar.bStart, hyperPar.BStart,
                                  hyperPar.cStart, hyperPar.CStart, par.mu, par.sigma,
                                  hyperPar.g, hyperPar.G, par.df, hyperPar.trans,
                                  hyperPar.a0, hyperPar.b0, hyperPar.d);
   if (K > 1)
   {
      /* Compute likelihood of Dirichlet prior */
      mixprior += priormixlik_dirichlet(par.weight,
                                        hyperPar.weightStart);
      cdpost  = mixlik + mixprior + dataC.postS;
      entropy = dataC.entropy;
   }
}
#endif /* __FINMIX_LOGSTUDENTIND_H__ */



