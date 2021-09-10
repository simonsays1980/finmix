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

#ifndef __FINMIX_PARSTUDENTFIX_H__
#define __FINMIX_PARSTUDENTFIX_H__

#include "FinmixModel.h"
#include "PriorStudentFix.h"

class ParStudentFix {
public:
arma::rowvec mu;
arma::rowvec sigma;
arma::rowvec df;
arma::rowvec acc;
bool INDEPENDENT;

ParStudentFix (const bool&, const FinmixModel&);
virtual ~ParStudentFix ()
{
}
virtual void update(const PriorStudentFix&);
virtual void permute(const arma::urowvec&,
                     const arma::urowvec&);
};
#endif /* __FINMIX_PARSTUDENTFIX_H__ */



