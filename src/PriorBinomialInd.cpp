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
#include "PriorBinomialInd.h"
#include "ParBinomialInd.h"
#include "posterior.h"

// =============================================================
// Constructor
// -------------------------------------------------------------

/**
 * -------------------------------------------------------------
 * PriorBinomialInd
 * @brief   Constructs a PriorBinomialInd object from a model.
 * @par prior   FinmixPrior object holding prior info
 * @return  PriorBinomialInd object
 * @detail  The only difference to a PriorBinomialFix object is
 *          the weight vector, all other members stay the same.
 *          This is achieved by a virtual inheritance from the
 *          PriorBinomialFix class.
 * @see FinmixPrior, PriorBinomialFix
 * @author Lars Simon Zehnder
 * ------------------------------------------------------------
 **/
PriorBinomialInd::PriorBinomialInd (const FinmixPrior& prior) :
   PriorBinomialFix(prior),
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
 * @par hyperPar    ParBinomialInd object containing the model
 *                  parameters
 * @detail  Updates the hyper parameters by computing posterior
 *          parameters for a Gamma prior for the component par-
 *          ameters and a Dirchlet prior for the weights.
 *          For updating the prior of the component parameters
 *          it is made use of the inheritance scheme and the
 *          corresponding update member function of the
 *          ParBinomialFix class is called.
 * @see PriorBinomialFix::update, ParBinomialInd,
 *      posterior_multinomial
 * @author  Lars Simon Zehnder
 * ------------------------------------------------------------
 **/
void PriorBinomialInd::update(const unsigned int& K, const arma::mat& y,
                              arma::ivec& S, const arma::vec& T, const ParBinomialInd& par)
{
   PriorBinomialFix::update(K, y, S, T, par);
   weightPost = posterior_multinomial(K, S, weightStart);
}

