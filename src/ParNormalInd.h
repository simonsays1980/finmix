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

#ifndef __FINMIX_PARNORMALIND_H__
#define __FINMIX_PARNORMALIND_H__

#include "ParNormalFix.h"
#include "PriorNormalInd.h"

class ParNormalInd : public ParNormalFix {
public:
arma::rowvec weight;

ParNormalInd (const bool&,
              const FinmixModel&);
virtual ~ParNormalInd ()
{
}
using ParNormalFix::update;
void update(const PriorNormalInd&);
};
#endif /* __FINMIX_PARNORMALIND_H__ */



