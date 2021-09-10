/******************************************************************************
*
* Copyright (C) 2013 Lars Simon Zehnder. All Rights Reserved.
*
* Author: Lars Simon Zehnder <simon.zehnder@gmail.com>
*
* This file is part of the R package 'finmix'.
*
* 'finmix' is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published
* by the Free Software Foundatio, either version 3 of the License, or
* any later version.
*
* 'finmix' is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with 'finmix'. If not, see <http://www.gnu.org/licenses/>.
*
******************************************************************************/


#ifndef ADAPTER_H
#define ADAPTER_H

#include <RcppArmadillo.h>
#include "BASE.h"

// =============================================================
// ADAPTER class (to be reviewed)
// -------------------------------------------------------------
/* @brief   Is used as the root for any layer combination.
 * @detail  This is the outer wrapper for a interlaced mixin
 *          layer construct. It defines a constructor for
 *          such that all necessary parameters can be provided.
 *          It inherits directly from the base class BASE and
 *          from any mxin layer above defined by 'Super'.
 *          Note, that the ADAPTER has actually no extra 'Node'
 *          and 'Output' mixin defined. It just takes already
 *          defined (or better refined) inner mixins 'Node' and
 *          'Output' from its Super class.
 * @see BASE, HIER, POST, IND, FIX
 * @author  Lars Simon Zehnder
 *
 * ============================================================
 * @review  An adapter class is in this setting probably not
 *          needed as all mixin layers have the same default
 *          parameters for their constructors respectively.
 *          Therefore any interlacing with no restruction in
 *          ordering can be done.
 * ------------------------------------------------------------
 **/
template <typename Super>
class ADAPTER : public Super, public BASE {
public:
ADAPTER ()
{
}
ADAPTER (const FinmixData&, const FinmixModel&, const
         FinmixPrior&, const FinmixMCMC&, Rcpp::S4&);
virtual void update();
virtual void store(const unsigned int&);
};

/**
 * ------------------------------------------------------------
 * ADPATER<Super>::ADAPTER
 * @brief   Constructs an ADAPTER object of any type 'Super'
 *          given as parameter to template.
 * @par data    a FinmixData object to hold all data
 * @par model   a FinmixModel object to hold all model
 *              information
 * @par prior   a FinmixPrior object to hold any information
 *              about the model prior
 * @par mcmc    a FinmixMCMC object to hold any configuration
 *              parameters for the Gibbs sampling algorithm
 * @par classS4     a Rcpp::S4 class object containing all
 *                  containers for output storage.
 * @detail  This is actually the main part of the ADAPTER. The
 *          constructor of the ADAPTER template contains all
 *          parameters needed to construct any upper mixin
 *          layers in an application. This constructor makes
 *          arbitrary interlacing of mixin layers possible.
 * @see FIX, HIER, IND, POST, BASE
 * @author Lars Simon Zehnder
 * -----------------------------------------------------------
 **/
template <typename Super>
ADAPTER <Super>::ADAPTER (const FinmixData& data, const FinmixModel& model, const
                          FinmixPrior& prior, const FinmixMCMC& mcmc, Rcpp::S4& classS4) :
   Super(data, model, prior, mcmc, classS4), BASE()
{
}

/**
 * -------------------------------------------------------
 * ADAPTER<Super>::update
 * @brief   Triggers the update process for each step
 *          the sampler. Passes responsibility to 'Node's
 *          'update()' method.
 * @see Super::Node::update
 * @author  Lars Simon Zehnder
 * -------------------------------------------------------
 **/
template <typename Super>
void ADAPTER <Super>::update()
{
   Super::update();
}

/**
 * -------------------------------------------------------
 * ADAPTER<Super>::store
 * @brief   Triggers the store process for each step of
 *          the sampler. Passes responsibility to 'Output's
 *          'store()' method.
 * @see Super::Output::store
 * @author  Lars Simon Zehnder
 * -------------------------------------------------------
 **/
template <typename Super>
void ADAPTER <Super>::store(const unsigned int& m)
{
   Super::store(m);
}
#endif
