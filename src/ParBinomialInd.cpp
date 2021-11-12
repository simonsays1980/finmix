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
#include "ParBinomialInd.h"

// =============================================================
// Constructor
// -------------------------------------------------------------

/**
 * -------------------------------------------------------------
 * ParBinomialInd
 * @brief   Constructs a ParBinomialInd object from a model.
 * @par STARTPAR    indicator for starting with sampling the
 *                  parameters
 * @par model       FinmixModel object holding model info
 * @return  ParBinomialInd object
 * @detail  The only difference to a ParBinomialFix object is
 *          the weight vector, all other members stay the same.
 *          This is achieved by a virtual inheritance from the
 *          ParBinomialFix class.
 * @see FinmixModel, ParBinomialFix
 * @author Lars Simon Zehnder
 * ------------------------------------------------------------
 **/
ParBinomialInd::ParBinomialInd (const bool& STARTPAR,
                                const FinmixModel& model) : ParBinomialFix(STARTPAR, model)
{
   if (!STARTPAR && model.K > 1)
   {
      weight = model.weight;
   }
}

// ============================================================
// Update
// ------------------------------------------------------------

/**
 * ------------------------------------------------------------
 * @brief   Updates the model parameters.
 * @par hyperPar    PriorBinomialInd object containing the hyper
 *                  parameters
 * @detail  Updates the parameters by sampling from a Beta
 *          distribution for the component parameters and a
 *          Dirichlet distribution for the weights. All hyper
 *          üparameters are provided by the PriorBinomialInd
 *          argument. For updating the component parameters it
 *          is made use of the inheritance scheme and the
 *          corresponding update member function of the
 *          ParBinomialFix class is called.
 * @see ParBinomialFix::update, PriorBinomialInd, rdirichlet
 * @author  Lars Simon Zehnder
 * ------------------------------------------------------------
 **/
void ParBinomialInd::update(const PriorBinomialInd& hyperPar) override
{
   ParBinomialFix::update(hyperPar);
   weight = rdirichlet(hyperPar.weightPost);
}
