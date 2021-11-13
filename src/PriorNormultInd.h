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

#ifndef __FINMIX_PRIORNORMULTIND_H__
#define __FINMIX_PRIORNORMULTIND_H__

#include "PriorNormultFix.h"

/* Forward declaration */
class ParNormultInd;

class PriorNormultInd : public PriorNormultFix {
public:
arma::rowvec weightStart;
arma::rowvec weightPost;

PriorNormultInd (const FinmixPrior&);
virtual ~PriorNormultInd ()
{
}
using PriorNormultFix::update;
void update(const unsigned int&,
                    const arma::mat&, arma::ivec&,
                    const arma::vec&, ParNormultInd&);
};
#endif /* __FINMIX_PRIORNORMULTIND_H__ */



