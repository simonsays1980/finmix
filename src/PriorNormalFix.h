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

#ifndef __FINMIX_PRIORNORMALFIX_H__
#define __FINMIX_PRIORNORMALFIX_H__

#include "FinmixPrior.h"

/* Forward declaration */
class ParNormalFix;
class PriorNormalFix {
    public:
        arma::rowvec bStart;
        arma::rowvec BStart;
        arma::rowvec cStart;
        arma::rowvec CStart;
        arma::rowvec bPost;
        arma::rowvec BPost;
        arma::rowvec cPost;
        arma::rowvec CPost;
        const bool HIER;
        const bool INDEPENDENT;
        double g;
        double G;

        PriorNormalFix ();
        PriorNormalFix (const FinmixPrior&);
        virtual ~PriorNormalFix () {}
        virtual void update (const unsigned int&,
                const arma::mat&, arma::ivec&, 
                const arma::vec&, ParNormalFix&);
        virtual void updateHier (const ParNormalFix&);        
};
#endif /* __FINMIX_PRIORNORMALFIX_H__ */



