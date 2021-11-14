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

#ifndef __FINMIX_LOGNORMALIND_H__
#define __FINMIX_LOGNORMALIND_H__

#include "LogNormalFix.h"
#include "ParNormalInd.h"
#include "PriorNormalInd.h"

class LogNormalInd : public LogNormalFix {
public:
double cdpost;
double entropy;
double maxcdpost;

LogNormalInd ();
virtual ~LogNormalInd ()
{
}
void update(const unsigned int&, const arma::mat&,
            arma::ivec&, const arma::mat&, const arma::vec&,
            const ParNormalInd&, const PriorNormalInd&);
};

LogNormalInd::LogNormalInd () : LogNormalFix(),
   cdpost(0.0), entropy(0.0), maxcdpost(0.0)
{
}

void LogNormalInd::update(const unsigned int& K,
                          const arma::mat& y, arma::ivec& S, const arma::mat& expos,
                          const arma::vec& T, const ParNormalInd& par,
                          const PriorNormalInd& hyperPar)
{
   liklist   lik   = likelihood_normal(y, par.mu, par.sigma);
   DataClass dataC = classification(S, lik, par.weight);

   S      = dataC.newS;
   mixlik = dataC.mixLik;
   /* Compute likelihood of mixture prior */
   mixprior = priormixlik_normal(hyperPar.INDEPENDENT,
                                 hyperPar.HIER, hyperPar.bStart, hyperPar.BStart,
                                 hyperPar.cStart, hyperPar.CStart, par.mu,
                                 par.sigma, hyperPar.g, hyperPar.G);
   if (K > 1)
   {
      /* Compute likelihood of Dirichlet prio */
      mixprior += priormixlik_dirichlet(par.weight,
                                        hyperPar.weightStart);
      cdpost  = mixlik + mixprior + dataC.postS;
      entropy = dataC.entropy;
   }
}
#endif /* __FINMIX_LOGNORMALIND_H__ */



