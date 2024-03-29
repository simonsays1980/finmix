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
#ifndef __FINMIX_PARCONDPOISSONIND_H_
#define __FINMIX_PARCONDPOISSONIND_H_

#include <RcppArmadillo.h>
#include "ParCondPoissonFix.h"
#include "PriorCondPoissonInd.h"

class ParCondPoissonInd : public ParCondPoissonFix {
public:
arma::rowvec weight;

ParCondPoissonInd (const bool&,
                   const FinmixModel&);
virtual ~ParCondPoissonInd ()
{
}
using ParCondPoissonFix::update;
void update(const PriorCondPoissonInd&);
};
#endif // __FINMIX_PARCONDPOISSONIND_H_
