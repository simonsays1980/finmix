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

#ifndef __FINMIX_PRIORSTUDMULTFIX_H__
#define __FINMIX_PRIORSTUDMULTFIX_H__

#include <RcppArmadillo.h>
#include "FinmixPrior.h"

/* Forward declaration */
class ParStudmultFix;

class PriorStudmultFix {
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

arma::rowvec mhTune;
const bool HIER;
bool INDEPENDENT;
std::string dftype;
double g;
arma::mat G;
double trans;
double a0;
double b0;
double d;

PriorStudmultFix ();
PriorStudmultFix (const FinmixPrior&);
~PriorStudmultFix ()
{
}
void update(const unsigned int&,
            const arma::mat&, arma::ivec&,
            const arma::vec&, ParStudmultFix&);
void updateDf(const unsigned int& K,
              const arma::mat& y, const arma::ivec& S,
              ParStudmultFix& par);
void updateHier(const ParStudmultFix&);
};
#endif /* __FINMIX_PRIORSTUDMULTFIX_H__ */



