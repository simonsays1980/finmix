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

#ifndef __FINMIX_PRIORSTUDENTFIX_H__
#define __FINMIX_PRIORSTUDENTFIX_H__

#include "FinmixPrior.h"

/* Forward declaration */
class ParStudentFix;
class PriorStudentFix {
    public:
        arma::rowvec bStart;
        arma::rowvec BStart;
        arma::rowvec cStart;
        arma::rowvec CStart;
        arma::rowvec bPost;
        arma::rowvec BPost;
        arma::rowvec cPost;
        arma::rowvec CPost;
        arma::rowvec mhTune;        
        const bool HIER;
        const bool INDEPENDENT;
        std::string dftype;
        double g;
        double G;
        double trans;
        double a0;
        double b0;
        double d;

        PriorStudentFix ();
        PriorStudentFix (const FinmixPrior&);
        virtual ~PriorStudentFix () {}
        virtual void update (const unsigned int&, 
                const arma::mat&, arma::ivec&,
                const arma::vec&, ParStudentFix&);
        virtual void updateDf (const unsigned int&, const arma::mat&,
                const arma::ivec&, ParStudentFix&);
        virtual void updateHier (const ParStudentFix&);
};
#endif /* __FINMIX_PRIORSTUDENTFIX_H__ */



