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

#ifndef __FINMIX_PRIORNORMALIND_H__
#define __FINMIX_PRIORNORMALIND_H__

#include "PriorNormalFix.h"

/* Forward declaration */
class ParNormalInd;

class PriorNormalInd : public PriorNormalFix {
public:
arma::rowvec weightStart;
arma::rowvec weightPost;

PriorNormalInd (const FinmixPrior&);
virtual ~PriorNormalInd ()
{
}
void update(const unsigned int&,
            const arma::mat&, arma::ivec&,
            const arma::vec&, ParNormalInd&);
};
#endif /* __FINMIX_PRIORNORMALIND_H__ */



