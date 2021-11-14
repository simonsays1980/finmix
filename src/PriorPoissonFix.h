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
* by the Free Software Foundatio, either version 3 of the License, or
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
#ifndef PRIORPOISSONFIX_H
#define PRIORPOISSONFIX_H

#include <RcppArmadillo.h>
#include "FinmixPrior.h"

class ParPoissonFix;
class PriorPoissonFix {
public:
arma::rowvec aStart;
arma::rowvec bStart;
arma::rowvec aPost;
arma::rowvec bPost;
const bool HIER;
double g;
double G;

PriorPoissonFix ();
PriorPoissonFix (const FinmixPrior&);
virtual ~PriorPoissonFix ()
{
}
virtual void update(const unsigned int&,
                    const arma::mat&, arma::ivec&,
                    const arma::vec&, const ParPoissonFix&);
virtual void updateHier(const ParPoissonFix&);
};
#endif
