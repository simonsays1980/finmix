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

#ifndef __FINMIX_PARNORMULTFIX_H__
#define __FINMIX_PARNORMULTFIX_H__

#include "FinmixModel.h"
#include "PriorNormultFix.h"

class ParNormultFix {
public:
arma::mat mu;
arma::cube sigma;
arma::cube sigmainv;
bool INDEPENDENT;

ParNormultFix (const bool&, const FinmixModel&);
virtual ~ParNormultFix ()
{
}
virtual void update(PriorNormultFix&);
virtual void permute(const arma::urowvec&,
                     const arma::urowvec&);
};
#endif /* __FINMIX_PARNORMULTFIX_H__ */



