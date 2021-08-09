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
 * by the Free Software Foundatio, either version 3 of the License, or
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
#include "ParCondPoissonInd.h"

ParCondPoissonInd::ParCondPoissonInd (const bool& STARTPAR, 
	const FinmixModel& model) :
		ParCondPoissonFix(STARTPAR, model), 
		weight(model.K) 
{
	if (!STARTPAR && model.K > 1) {
		weight = model.weight;
	}
}

void ParCondPoissonInd::update (const PriorCondPoissonInd& hyperPar)
{
	ParCondPoissonFix::update(hyperPar);
	weight = rdirichlet(hyperPar.weightPost);
}
