/******************************************************************************
 *
 * Copyright (C) 2013 Lars Simon Zehnder. All Rights Reserved.
 *
 * Author: Lars Simon Zehnder <simon.zehnder@gmail.com>
 *
 * This file is part of the R package finmix.
 *
 * finmix is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * any later version.
 *
 * finmix is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with finmix. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#ifndef __FINMIX_PAREXPONENTIALFIX_H__
#define __FINMIX_PAREXPONENTIALFIX_H__

#include "FinmixModel.h"
#include "PriorExponentialFix.h"
#include "distributions.h"

class ParExponentialFix {
	public: 
		arma::rowvec lambda;
		
		ParExponentialFix (const bool&, 
				const FinmixModel&);
		virtual ~ParExponentialFix () {}
		void update (const PriorExponentialFix&);
        virtual void permute (const arma::urowvec&, 
                const arma::urowvec&);
};
#endif
