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
#include "ParCondPoissonFix.h"
#include "rtruncnorm.h"
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
 * @return          an object of class ParPoissonFix
 * @detail  If STARTPAR == FALSE it should be started by sampling
 *          the indicators and starting parameters are provided 
 *          by the model parameter
 * @see ?model in R
 * @author Lars Simon Zehnder
 * -------------------------------------------------------------
 **/

ParCondPoissonFix::ParCondPoissonFix (const bool& STARTPAR, 
		const FinmixModel& model) : lambda(model.K),
    acc(0.0)
{
	if(!STARTPAR && model.K > 1) {
		arma::rowvec tmp = Rcpp::as<arma::rowvec>
				((SEXP) model.par["lambda"]);
		lambda = tmp;
	}
} 

double ParCondPoissonFix::metropolis (const arma::rowvec& can, 
        const PriorCondPoissonFix& hyperPar) 
{
    unsigned int K = can.n_elem;
    double output = std::pow(can(0), hyperPar.Q(0)) 
        * std::exp(-hyperPar.N(0) * can(0))
        * R::dunif(can(0), hyperPar.a, hyperPar.b, 0);
    for(unsigned int k = 1; k < K; ++k) {
        output *= std::pow(can(k), hyperPar.Q(k)) 
            * std::exp(-hyperPar.N(k) * can(k))
            * do_dtruncnorm(can(k), can(k - 1), R_PosInf, can(k - 1), hyperPar.s);        
    }
    return output;
}

double ParCondPoissonFix::proposal_prob (const arma::rowvec& can,
        const PriorCondPoissonFix& hyperPar) 
{
    const unsigned int K = can.n_elem;
    double proposal = R::dunif(can(0), hyperPar.a, hyperPar.b, 0);    
    for(unsigned int k = 1; k < K; ++k) {
        proposal *= do_dtruncnorm(can(k), can(k - 1), R_PosInf, can(k - 1), hyperPar.s);
    }
    return proposal;
}
// =============================================================
// Update
// -------------------------------------------------------------

/** 
 * -------------------------------------------------------------
 * update
 * @brief   Updates the parameters of the conditional Poisson model
 * @par hyperPar    object of class PriorCondPoissonFix, holds
 *                  hyper parameters for sampling.
 * @details draws samples from a Gamma prior
 * @see ?prior in R
 * @author Lars Simon Zehnder
 * -------------------------------------------------------------
 **/
void ParCondPoissonFix::update (const PriorCondPoissonFix& hyperPar) 
{
    acc = 0.0;
    unsigned int K = lambda.n_elem;
    arma::rowvec can(K);
    can(0) = R::runif(hyperPar.a, hyperPar.b);
    for(unsigned int k = 1; k < K; ++k) {
        can(k) = do_rtruncnorm(1, can(k - 1), R_PosInf, can(k - 1), hyperPar.s)(0);
    }
    double aprob = std::min(1.0, (metropolis(can, hyperPar) / metropolis(lambda, hyperPar))
            / (proposal_prob(can, hyperPar) / proposal_prob(lambda, hyperPar)));

/*    double aprob = std::min(1.0, (metropolis(can, hyperPar.N, hyperPar.Q, hyperPar.s, hyperPar.a, hyperPar.b)
                /metropolis(lambda, hyperPar.N, hyperPar.Q, hyperPar.s, hyperPar.a, hyperPar.b))
            / ((R::dunif(can(1), hyperPar.a, hyperPar.b, 0) * do_dtruncnorm(can(0), can(1), R_PosInf, can(1), hyperPar.s))
            / (R::dunif(lambda(1), hyperPar.a, hyperPar.b, 0) * do_dtruncnorm(lambda(0), can(1), R_PosInf, can(1), hyperPar.s))));*/
    double u = R::runif(0.0, 1.0);
    if(u < aprob && NA_REAL != arma::prod(can)) {        
        lambda = can;
        acc = 1.0; 
    }
}

void ParCondPoissonFix::permute (const arma::urowvec& compIndex,
        const arma::urowvec& permIndex) {}
