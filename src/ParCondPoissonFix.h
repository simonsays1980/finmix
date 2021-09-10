/******************************************************************************
*
* Copyright (C) 2013 Lars Simon Zehnder. All Rights Reserved.
*
* Author: Lars Simon Zehnder <simon.zehnder@gmail.com>
*
* This file is part of the R package finmix.
*
* finmix is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published
* by the Free Software Foundation, either version 3 of the License, or
* any later version.
*
* finmix is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with finmix. If not, see <http://www.gnu.org/licenses/>.
*
******************************************************************************/
#ifndef __FINMIX_PARCONDPOISSONFIX_H_
#define __FINMIX_PARCONDPOISSONFIX_H_

#include <RcppArmadillo.h>
#include "FinmixModel.h"
#include "PriorCondPoissonFix.h"
#include "distributions.h"

class ParCondPoissonFix {
public:
arma::rowvec lambda;
double acc;

ParCondPoissonFix (const bool&,
                   const FinmixModel&);
virtual ~ParCondPoissonFix ()
{
}
double metropolis(const arma::rowvec&,
                  const PriorCondPoissonFix&);
double proposal_prob(const arma::rowvec&,
                     const PriorCondPoissonFix&);
void update(const PriorCondPoissonFix&);
virtual void permute(const arma::urowvec&,
                     const arma::urowvec&);
};
#endif // __FINMIX_PARCONDPOISSON_H_
