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

#ifndef __FINMIX_PARNORMALFIX_H__
#define __FINMIX_PARNORMALFIX_H__

#include "FinmixModel.h"
#include "PriorNormalFix.h"

class ParNormalFix {
public:
arma::rowvec mu;
arma::rowvec sigma;
bool INDEPENDENT;

ParNormalFix (const bool&, const FinmixModel&);
virtual ~ParNormalFix ()
{
}
void update(const PriorNormalFix&);
virtual void permute(const arma::urowvec&,
                     const arma::urowvec&);
};
#endif /* __FINMIX_PARNORMALFIX_H__ */



