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

#ifndef __FINMIX_PARSTUDENTIND_H__
#define __FINMIX_PARSTUDENTIND_H__

#include "ParStudentFix.h"
#include "PriorStudentInd.h"

class ParStudentInd : virtual public ParStudentFix {
    public:
        arma::rowvec weight;

        ParStudentInd (const bool&, 
                const FinmixModel&);
        virtual ~ParStudentInd () {}
        virtual void update (const PriorStudentInd&);
};
#endif /* __FINMIX_PARSTUDENTIND_H__ */



