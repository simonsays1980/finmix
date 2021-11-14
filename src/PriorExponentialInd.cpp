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
#include "PriorExponentialInd.h"
#include "ParExponentialInd.h"
#include "posterior.h"

// =============================================================
// Constructor
// -------------------------------------------------------------

/**
 * -------------------------------------------------------------
 * PriorExponentialInd
 * @brief   Constructs a PriorExponentialInd object from a model.
 * @par prior   FinmixPrior object holding prior info
 * @return  PriorExponentialInd object
 * @detail  The only difference to a PriorExponentialFix object is
 *          the weight vector, all other members stay the same.
 *          This is achieved by a virtual inheritance from the
 *          PriorExponentialFix class.
 * @see FinmixPrior, PriorExponentialFix
 * @author Lars Simon Zehnder
 * ------------------------------------------------------------
 **/
PriorExponentialInd::PriorExponentialInd (const FinmixPrior& prior) :
   PriorExponentialFix(prior),
   weightStart(prior.weight),
   weightPost(prior.weight)
{
}

// ============================================================
// Update
// ------------------------------------------------------------

/**
 * ------------------------------------------------------------
 * @brief   Updates the hyper parameters.
 * @par hyperPar    ParExponentialInd object containing the model
 *                  parameters
 * @detail  Updates the hyper parameters by computing posterior
 *          parameters for a Gamma prior for the component par-
 *          ameters and a Dirchlet prior for the weights.
 *          For updating the prior of the component parameters
 *          it is made use of the inheritance scheme and the
 *          corresponding update member function of the
 *          ParExponentialFix class is called.
 * @see PriorExponentialFix::update, ParExponentialInd,
 *      posterior_multinomial
 * @author  Lars Simon Zehnder
 * ------------------------------------------------------------
 **/
void PriorExponentialInd::update(const unsigned int& K, const arma::mat& y,
                                 arma::ivec& S, const arma::vec& T, const ParExponentialInd& par)
{
   PriorExponentialFix::update(K, y, S, T, par);
   weightPost = posterior_multinomial(K, S, weightStart);
}
