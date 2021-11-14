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

#ifndef __FINMIX_PARNORMULTIND_H__
#define __FINMIX_PARNORMULTIND_H__

#include "ParNormultFix.h"
#include "PriorNormultInd.h"

class ParNormultInd : public ParNormultFix {
public:
arma::rowvec weight;

ParNormultInd (const bool&, const FinmixModel&);
~ParNormultInd ()
{
}
using ParNormultFix::update;
void update(PriorNormultInd&);
};
#endif /* __FINMIX_PARNORMULTIND_H__ */



