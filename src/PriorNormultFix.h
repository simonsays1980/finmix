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

#ifndef __FINMIX_PRIORNORMULTFIX_H__
#define __FINMIX_PRIORNORMULTFIX_H__

#include <RcppArmadillo.h>
#include "FinmixPrior.h"

/* Forward declaration */
class ParNormultFix;

class PriorNormultFix {
public:
arma::mat bStart;
arma::cube BStart;
arma::cube BInvStart;
arma::rowvec N0Start;
arma::rowvec cStart;
arma::cube CStart;

arma::mat bPost;
arma::cube BPost;
arma::cube BInvPost;
arma::rowvec N0Post;
arma::rowvec cPost;
arma::cube CPost;
arma::rowvec logdetC;
const bool HIER;
bool INDEPENDENT;
double g;
arma::mat G;

PriorNormultFix ();
PriorNormultFix (const FinmixPrior&);
~PriorNormultFix ()
{
}
void update(const unsigned int&,
            const arma::mat&, arma::ivec&,
            const arma::vec&, ParNormultFix&);
void updateHier(const ParNormultFix&);
};
#endif /* __FINMIX_PRIORNORMULTFIX_H__ */



