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

#ifndef __FINMIX_PRIORBINOMIALFIX_H__
#define __FINMIX_PRIORBINOMIALFIX_H__

#include "FinmixPrior.h"

/* Forward declaration */
class ParBinomialFix;

class PriorBinomialFix {
public:
arma::rowvec aStart;
arma::rowvec bStart;
arma::rowvec aPost;
arma::rowvec bPost;

PriorBinomialFix ()
{
}
PriorBinomialFix (const FinmixPrior&);
virtual ~PriorBinomialFix ()
{
}
virtual void update(const unsigned int&,
                    const arma::mat&, arma::ivec&,
                    const arma::vec&, const ParBinomialFix&);
/* Needs to be defined even if not used
 * as FIX::update() calls it regularly
 */
virtual void updateHier(const ParBinomialFix&)
{
};
};
#endif /* __FINMIX_PRIORBINOMIALFIX_H__ */



