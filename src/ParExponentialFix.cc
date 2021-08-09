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
 * along with 'finmix'. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#include "ParExponentialFix.h"

// =============================================================
// Constructor
// -------------------------------------------------------------

/**
 * -------------------------------------------------------------
 * @brief   Constructs object from model parameters.
 * @par STARTPAR    boolean, indicating if it should be started
 *                  by sampling the parameters
 * @par model       FinmixModel object, holding model 
 *                  definitions and starting parameters
 * @return          an object of class ParExponentialFix
 * @detail  If STARTPAR == FALSE it should be started by sampling
 *          the indicators and starting parameters are provided 
 *          by the model parameter
 * @see ?model in R
 * @author Lars Simon Zehnder
 * -------------------------------------------------------------
 **/

ParExponentialFix::ParExponentialFix (const bool& STARTPAR, 
		const FinmixModel& model) : lambda(model.K) 
{
	if(!STARTPAR && model.K > 1) {
		arma::rowvec tmp = Rcpp::as<arma::rowvec>
				((SEXP) model.par["lambda"]);
		lambda = tmp;
	}
} 

// =============================================================
// Update
// -------------------------------------------------------------

/** 
 * -------------------------------------------------------------
 * update
 * @brief   Updates the parameters of the Exponential model
 * @par hyperPar    object of class PriorExponentialFix, holds
 *                  hyper parameters for sampling.
 * @details draws samples from a Gamma prior
 * @see ?prior in R
 * @author Lars Simon Zehnder
 * -------------------------------------------------------------
 **/

void ParExponentialFix::update (const PriorExponentialFix& hyperPar) 
{
	lambda = rgammaprod(hyperPar.aPost, hyperPar.bPost);
}

// =============================================================
// Permute
// -------------------------------------------------------------

/**
 * -------------------------------------------------------------
 * permute
 * @brief   Permutes the parameters of the Exponential model
 * @param   compIndex   row vector with indices 1 to K
 * @param   permIndex   row vector with permuted indices
 * @see ?mcmc
 * @author  Lars Simon Zehnder
 * -------------------------------------------------------------
 **/
void ParExponentialFix::permute (const arma::urowvec& compIndex,
        const arma::urowvec& permIndex)
{
    lambda(compIndex) = lambda(permIndex);
}
