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
#include "ParBinomialFix.h"

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
 * @return          an object of class ParBinomialFix
 * @detail  If STARTPAR == FALSE it should be started by sampling
 *          the indicators and starting parameters are provided 
 *          by the model parameter
 * @see ?model in R
 * @author Lars Simon Zehnder
 * -------------------------------------------------------------
 **/
ParBinomialFix::ParBinomialFix(const bool& STARTPAR,
        const FinmixModel& model) : p(model.K)
{
    if (!STARTPAR && model.K > 1) {
        arma::rowvec tmp = Rcpp::as<arma::rowvec>
            ((SEXP) model.par["p"]);
        p = tmp;
    }
}

// =============================================================
// Update
// -------------------------------------------------------------

/** 
 * -------------------------------------------------------------
 * update
 * @brief   Updates the parameters of the Binomial model
 * @par hyperPar    object of class PriorBinomialFix, holds
 *                  hyper parameters for sampling.
 * @details draws samples from a Beta prior
 * @see ?prior in R
 * @author Lars Simon Zehnder
 * -------------------------------------------------------------
 **/
void ParBinomialFix::update (const PriorBinomialFix& hyperPar) 
{
    p = rbetaprod(hyperPar.aPost, hyperPar.bPost);
}

void ParBinomialFix::permute (const arma::urowvec& compIndex,
        const arma::urowvec& permIndex) 
{
    p(compIndex) = p(permIndex);
}
