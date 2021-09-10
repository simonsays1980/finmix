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
#ifndef __FINMIX_PRIOREXPONENTIALIND_H__
#define __FINMIX_PRIOREXPONENTIALIND_H__

#include "PriorExponentialFix.h"

class ParExponentialInd;
class PriorExponentialInd : virtual public PriorExponentialFix {
public:
arma::rowvec weightStart;
arma::rowvec weightPost;

PriorExponentialInd (const FinmixPrior&);
virtual ~PriorExponentialInd ()
{
}
virtual void update(const unsigned int&,
                    const arma::mat&, arma::ivec&,
                    const arma::vec&, const ParExponentialInd&);
};
#endif /* __FINMIX_PRIOREXPONENTIALIND_H__ */
