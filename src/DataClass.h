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

#ifndef DATACLASS_H
#define DATACLASS_H

#include <RcppArmadillo.h>
#include "likelihood.h" // methods for calculating likelihoods

/* structure to contain:
 * logpy    : log likelihood values for 1 to K
 * prob     : classification probability matrix
 * newS     : sampled classifications from classification probability matrix
 * loglikcd : posterior complete data likelihood, knowing classifications
 * mixlik   : posterior mixture likelihood
 * entropy  : entropy of the posterior classification probability distribution
 * postS    : posterior of sampled classifications
 */
struct DataClass
{
   arma::mat  logPy;
   arma::mat  prob;
   arma::ivec newS;
   arma::vec  logLikCd;
   double     mixLik;
   double     entropy;
   double     postS;
};

/**
 * This function is used for all mixtures
 *
 */
inline DataClass
classification(const arma::ivec &S, const liklist &lik,
               const arma::rowvec& weight)
{
   const unsigned int N     = S.n_elem;
   const unsigned int K     = weight.n_elem;
   double             postS = 0.0;
   DataClass          dataC = DataClass();
//    lik.lh.print("lh: ");
//    lik.llh.print("llh: ");
   /* if indicators are not fix, they are simulated */
   arma::mat  p_m(N, K);
   arma::ivec newS(N);

   /* only multinomial indicator model implemented */
   for (unsigned int i = 0; i < N; ++i)
   {
      p_m.row(i) = lik.lh.row(i) % weight;
   }
   arma::vec sump_v = arma::sum(p_m, 1);                                                // N x 1 matrix
   arma::vec lsump  = arma::log(sump_v) + lik.maxl;                                     // N x 1 matrix
   double    mixlik = arma::sum(lsump);                                                 // mixture likelihood

   p_m.each_col() /= sump_v;                                                            // classification probability matrix
                                                                                        // N x K matrix

   dataC.prob = p_m;
   /* simulate only if true mixture */
   if (K > 1)
   {
      /* simulate classifications from probability matrix p */
      arma::vec rnd(N);
      GetRNGstate();
      for (unsigned int i = 0; i < N; ++i)
      {
         rnd(i) = R::runif(0.0, 1.0);
      }
      PutRNGstate();
      arma::mat  rndM  = arma::repmat(rnd, 1, K);                                       // N x K matrix
      arma::mat  cumSP = arma::cumsum(p_m, 1);                                          // cumulate along rows, N x K matrix
      arma::umat ind   = (cumSP > rndM);
      rndM = arma::conv_to<arma::mat>::from(ind);                                       // logical N x K matrix
      newS = arma::conv_to<arma::ivec>::from(arma::sum(rndM, 1));                       // new classifications

      /* compute posterior log likelihood of S */
      arma::imat Sm    = arma::repmat(newS, 1, K);                                      // N x K matrix of S
      arma::imat compM = arma::ones<arma::imat>(N, K);
      for (unsigned int k = 0; k < K; ++k)
      {
         compM.col(k) = compM.col(k) * (k + 1);
      }
      ind = (Sm == compM);                                                              // logical N x K matrix
      arma::mat indDouble = arma::conv_to<arma::mat>::from(ind);
      arma::vec postSm    = arma::sum(indDouble % p_m, 1);                              // sum along rows
      postSm = arma::log(postSm);
      postS  = arma::sum(postSm);
   }


   /* calculate entropy */
   arma::mat  logp(N, K);
   arma::uvec col_index(1);

   for (unsigned int k = 0; k < K; ++k)
   {
      col_index(0) = k;
      arma::uvec zero_index = arma::find(p_m.col(k) == 0);
      arma::uvec index      = arma::find(p_m.col(k));
      logp.submat(zero_index, col_index).fill(-99.00);
      logp.submat(index, col_index) = arma::log(p_m.submat(index, col_index));
   }
   double entropy = (-1.0) * arma::accu(logp % p_m);

   dataC.logPy   = lik.llh;
   dataC.prob    = p_m;
   dataC.newS    = newS;
   dataC.mixLik  = mixlik;
   dataC.entropy = entropy;
   dataC.postS   = postS;

   return dataC;
}


/**
 * -------------------------------------------------------------
 * classification_fix
 * @brief   Computes the complete data mixture log- likelihood
 *          for a model with fixed indicators.
 * @par K       number of components
 * @par S       fixed indicators
 * @par liklist a liklist object holding the likelihoods and log
 *              log-likelihoods for each observation and com-
 *              ponent as well as the maximum likelihood over
 *              components.
 * @details computes the complete data likelihood by summing the
 *          log-likelihoods over the component of each ob-
 *          servation indicated by S.
 * @see liklist
 * @author Lars Simon Zehnder
 * -------------------------------------------------------------
 **/
inline DataClass
classification_fix(const unsigned int K, const arma::ivec& S,
                   const liklist& lik)
{
   arma::vec  loglikcd(K);
   arma::uvec col_index(1);

   if (K > 1)
   {
      for (unsigned int k = 0; k < K; ++k)
      {
         col_index(0) = k;
         arma::uvec index = arma::find(S == k);
         loglikcd(k) = arma::accu(lik.llh.submat(index, col_index));
      }
   }
   else        /* no true mixture */
   {
      arma::vec sump_v = sum(lik.lh, 1);
      arma::vec lsump  = arma::log(sump_v) + lik.maxl;
      double    mixlik = arma::sum(lsump);
      //TODO: check if better to fill all K entries with mixlik
      loglikcd(0, 0) = mixlik;
   }
   DataClass dataC;

   dataC.logPy    = lik.llh;
   dataC.logLikCd = loglikcd;
   return dataC;
}
#endif
