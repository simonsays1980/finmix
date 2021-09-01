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
#ifndef __FINMIX_PRIORLIKELIHOOD_H__
#define __FINMIX_PRIORLIKELIHOOD_H__

#include <RcppArmadillo.h>
#include <R.h> 			// to interface with R
#include <Rmath.h> 		// for using internal R C-functions
#include "likelihood.h"
#include "distributions.h"
#include "rtruncnorm.h"
#include "PriorCondPoissonFix.h"

/**
 * Evaluates the prior likelihood for the 
 * weights. This function is used in every
 * Gibbs sampling. 
 * 
 */

// =============================================================
// Dirichlet prior
// -------------------------------------------------------------

/**
 * -------------------------------------------------------------
 * @brief   Computes the mixture log-likelihood for a 
 *          Dirichlet distribution. 
 * @par eta         weight parameters
 * @par prior_par   hyper parameters of the Dirichlet 
 *                  distribution
 * @return  Prior mixture log-likelihood
 * @detail  The log-lieklihood is computed via the density
 *          of the Dirichlet distribution evaluated at the 
 *          weights with hyper parameters determined in 
 *          prior_par
 * @author Lars Simon Zehnder
 * -------------------------------------------------------------
 **/
inline double 
priormixlik_dirichlet(const arma::rowvec &eta, 
			const arma::rowvec &prior_par) 
{
	unsigned int K      = eta.n_elem;
	double priormixlik  = 0.0;	
	/* Evaluate Dirichlet loglik */
	priormixlik = R::lgammafn(arma::accu(prior_par));
	for(unsigned int k = 0; k < K; ++k) {
		priormixlik += (prior_par(k) - 1) * std::log(eta(k));
        priormixlik -= R::lgammafn(prior_par(k));		 
	}
	return priormixlik;
}
// =============================================================
// Poisson prior
// -------------------------------------------------------------

/**
 * -------------------------------------------------------------
 * @brief   Computes the prior mixture log-likelihood for a 
 *          Poisson distribution. 
 * @par lambda      Poisson parameters, 1 x K 
 * @par prior_parA  Gamma hyper shape parameters, 1 x K
 * @par prior_parB  Gamma hyper rate parameters, 1 x K
 * @par hier        boolean value indicating if a hierarchical
 *                  prior is used
 * @par g           Gamma hyper shape parameter of the hierar-
 *                  chical prior
 * @par G           Gamma hyper rate parameter of the hierar-
 *                  chical prior
 * @return  the prior mixture log-likelihood value
 * @author Lars Simon Zehnder
 * ------------------------------------------------------------
 **/
inline double
priormixlik_poisson(const arma::rowvec& lambda, 
			const arma::rowvec& prior_parA,
			const arma::rowvec& prior_parB, 
			const bool &hier, const double &g, 
			const double &G) {

	unsigned int K = lambda.n_elem;
	double priormixlik = 0.0;
	if(!hier) {
		priormixlik = likelihood_gamma(lambda, prior_parA(0), prior_parB(0));
	}
	else { /* hierarchical prior */
		double gN = g + K * prior_parA(0); // prior_parA must be the start value.
		double GN = G + arma::accu(lambda);
		double b = gN/GN;
		double scale = 1.0/b;
		/* step 1: log likelihood of prior */
	 	for(unsigned int k = 0; k < K; ++k) {
			priormixlik += R::dgamma(lambda(k), prior_parA(k), scale, 1);
		}
			/**
		 * step 2: log likelihood of hyperprior with start 
		 * values of hyper parameters 
		 */
		scale = 1.0/G;
		priormixlik += R::dgamma(b, g, scale, 1);
	
		/**
		 * step 3: log likelihood of hyperprior with updated 
		 * hyper parameters 
 		 */
		scale = 1.0/GN;
		priormixlik -= R::dgamma(b, gN, scale, 1);

	}
	return priormixlik;	
}

inline double 
priormixlik_condpoisson (const arma::rowvec& lambda,
	const PriorCondPoissonFix& hyperPar)
{	
    unsigned int K = lambda.n_elem;
    double priormixlik = R::dunif(lambda(0), hyperPar.a, hyperPar.b, 1); 
    for(unsigned int k = 1; k < K; ++k) {
        priormixlik += std::log(do_dtruncnorm(lambda(k), lambda(k - 1), 
                    R_PosInf, lambda(k - 1), hyperPar.s));
    }
	return priormixlik;
}

// =============================================================
// Binomial distribution
// -------------------------------------------------------------

/**
 * -------------------------------------------------------------
 * @brief   Computes the prior mixture log-lieklihood for a 
 *          Binomial distribution. 
 * @par lambda      Poisson parameters, 1 x K 
 * @par prior_parA  Gamma hyper shape parameters, 1 x K
 * @par prior_parB  Gamma hyper rate parameters, 1 x K
 * @par hier        boolean value indicating if a hierarchical
 *                  prior is used
 * @par g           Gamma hyper shape parameter of the hierar-
 *                  chical prior
 * @par G           Gamma hyper rate parameter of the hierar-
 *                  chical prior
 * @return  the prior mixture log-likelihood value
 * @author Lars Simon Zehnder
 * ------------------------------------------------------------
 **/
inline
double priormixlik_binomial (const arma::rowvec& p,
        const arma::rowvec& prior_parA, 
        const arma::rowvec& prior_parB)
{
     const unsigned int K = p.n_elem;
     double priormixlik = 0.0;
     for (unsigned int k = 0; k < K; ++k) {
         priormixlik += (prior_parA(k) - 1.0) * std::log(p(k)) 
             + (prior_parB(k) - 1.0) * std::log(p(k));
         priormixlik -= R::lbeta(prior_parA(k), prior_parB(k));
     }
     return priormixlik;
}

inline
double priormixlik_normal (const bool& INDEPENDENT, const bool& HIER, 
        const arma::rowvec& bStart, const arma::rowvec& BStart,
        const arma::rowvec& cStart, const arma::rowvec& CStart,
        const arma::rowvec& mu, const arma::rowvec& sigma,
        const double& g, const double& G)
{
    const unsigned int K = mu.n_elem;
    double mixlik = 0.0 ;
    if ( INDEPENDENT ) {
        mixlik = arma::sum(arma::log(1.0 / BStart * M_PI * 2.0));
        mixlik += arma::sum(arma::pow(mu - bStart, 2.0) / (1.0 / BStart));
        mixlik *= -0.5;        
    } else { /* conditionally conjugate prior */
             /* here, B == N0 */
        mixlik = arma::sum(arma::log(sigma / BStart * M_PI * 2.0));
        mixlik += arma::sum(arma::pow(mu - bStart, 2.0) / (sigma / BStart));
        mixlik *= -0.5;
    }
    /* add likelihood for sigma */
    for (unsigned int k = 0; k < K; ++k) {
        mixlik += cStart(k) * std::log(CStart(k));
        mixlik -= R::lgammafn(cStart(k));
        mixlik -= CStart(k) / sigma(k);
        mixlik -= (cStart(k) + 1) * std::log(sigma(k));
    }
    if (HIER) {
        double gN = g + arma::sum(cStart);
        double GN = G + arma::sum( 1.0 / sigma );
        double Cstar = gN / GN;
        mixlik += R::dgamma(Cstar, g, 1.0 / G, 1);
        mixlik -= R::dgamma(Cstar, gN, 1.0 / GN, 1);
    }
    return mixlik;
}

inline
double priormixlik_normult (const bool& INDEPENDENT, const bool& HIER,
        const arma::mat& bStart, const arma::cube& BInvStart,const  arma::cube& BStart, 
        const arma::rowvec& cStart, arma::cube CStart, const arma::rowvec& logdetC,
        const double& g, const arma::mat& G, const arma::mat& mu, 
        const arma::cube& sigma)
{
    const unsigned int K = mu.n_cols;
    double mixlik = 0.0;
    if (INDEPENDENT) {
        mixlik += logdnormult(mu, bStart, BStart, BInvStart);
    } else { /* conditionally conjugate prior */
        mixlik += logdnormult(mu, bStart, BStart, BInvStart);
    }
    if (HIER) {
        arma::rowvec gvec(K);
        gvec.fill(g);
        double gN       = g  + arma::sum(cStart);
        arma::rowvec gNvec(K);
        gNvec.fill(gN);
        arma::mat GN    = G;
        for (unsigned int k = 0; k < K; ++k) {
            GN += CStart.slice(k);            
        }
        arma::mat Cstar = gN * arma::inv(GN);
        for (unsigned int k = 0; k < K; ++k) {
            CStart.slice(k) = Cstar;
        }
        mixlik += logdwishart(CStart, gvec, G, std::log(arma::det(G)));
        mixlik += logdwishart(CStart, gNvec, GN, std::log(arma::det(GN)));
    }
    /* Prior for sigma (Wishart for sigmainv) */
    mixlik += logdwishart(sigma, cStart, CStart, logdetC);    
    return mixlik;
}

inline
double priormixlik_student (const bool& INDEPENDENT, const bool& HIER,
        const arma::rowvec& bStart, const arma::rowvec BStart, 
        const arma::rowvec& cStart, arma::rowvec CStart,
        const arma::rowvec& mu, const arma::rowvec& sigma, 
        const double& g, const double& G,
        const arma::rowvec& df, const double& trans,
        const double& a0, const double& b0, const double& d)
{
    double loglik = priormixlik_normal(INDEPENDENT, HIER, bStart, BStart,
            cStart, CStart, mu, sigma, g, G);
    arma::rowvec fnu(mu.n_elem);
    arma::rowvec nu = df - trans;
    if (b0 == 1.0) {
        fnu = std::log(d + a0) + (a0 - 1) * arma::log(nu) 
            - (a0 + 1) * arma::log(nu + d); 
    } else {
        fnu = b0 * std::log(d) + (a0 - 1) * arma::log(nu)
            - (a0 + b0) * arma::log(nu + d) - R::lbeta(a0, b0);
    }
    loglik += arma::sum(fnu);
    return loglik;
}

inline 
double priormixlik_studmult (const bool& INDEPENDENT, const bool& HIER,
        const arma::mat& bStart, const arma::cube BInvStart, 
        const arma::cube& BStart, const arma::rowvec& cStart, 
        arma::cube CStart, const arma::rowvec& logdetC, 
        const double& g, const arma::mat& G, const arma::mat& mu,
        const arma::cube& sigma, const arma::rowvec& df, 
        const double& trans, const double& a0, const double& b0,
        const double& d)
{
    double loglik = priormixlik_normult(INDEPENDENT, HIER,
            bStart, BInvStart, BStart, cStart, CStart, logdetC, 
            g, G, mu, sigma); 
    arma::rowvec fnu(mu.n_elem);
    arma::rowvec nu = df - trans;
    if (b0 == 1.0) {
        fnu = std::log(d + a0) + (a0 - 1.0) * arma::log(nu)
            - (a0 + 1.0) * arma::log(nu + d);        
    } else {
        fnu = b0 * std::log(d) + (a0 - 1.0) * arma::log(nu)
            - (a0 + b0) * arma::log(nu + d) - R::lbeta(a0, b0);
    }   
    loglik += arma::sum(fnu);
    return loglik;
}
#endif
