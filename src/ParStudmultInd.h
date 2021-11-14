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

#ifndef __FINMIX_PARSTUDMULTIND_H__
#define __FINMIX_PARSTUDMULTIND_H__

#include "ParStudmultFix.h"
#include "PriorStudmultInd.h"

class ParStudmultInd : public ParStudmultFix {
public:
arma::rowvec weight;

ParStudmultInd (const bool&, const FinmixModel&);
~ParStudmultInd ()
{
}
using ParStudmultFix::update;
void update(PriorStudmultInd&);
};
#endif /* __FINMIX_PARSTUDMULTIND_H__ */



