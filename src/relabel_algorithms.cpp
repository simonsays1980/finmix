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

// [[Rcpp::depends(RcppArmadillo)]]

#include "algorithms.h"
#include "optimize.h"
#include "hungarian.h"
#include <nloptrAPI.h>
#include <cmath>


// ============================================================
// Stephens Relabeling Algorithm (1997a)
// ------------------------------------------------------------

/**
 * ------------------------------------------------------------
 * stephens1997a_poisson_cc
 * @brief   Defines Stephens (1997a) relabelling algorithm for
 *          Poisson models. The nlopt library is used for
 *          optimization (Nelder-Mead algorithm)
 * @par values1 sampled lambda parameters; M x K
 * @par values2 sampled weight parameters; M x K
 * @par pars    Gamma and Dirichlet hyper parameters
 * @par perm    matrix with all possible permutations of labels;
 * @return  matrix indicating the optimal labeling of sampled
 *          parameters; M x K
 * @detail  See Stephens (1997a)
 * @see nlopt
 * @author  Lars Simon Zehnder
 * ------------------------------------------------------------
 **/

// [[Rcpp::export]]

arma::imat stephens1997a_poisson_cc(Rcpp::NumericMatrix values1,
                                    Rcpp::NumericMatrix values2,
                                    arma::vec pars, const arma::umat perm)
{
   const unsigned int M          = values1.rows();
   const unsigned int K          = values1.cols();
   const unsigned int P          = perm.n_rows;
   const unsigned int n          = pars.n_elem;
   double             value      = 0.0;
   double             value_next = -10.0e-8;
   arma::mat          lambda(values1.begin(), M, K, true, true);
   arma::mat          weight(values2.begin(), M, K, true, true);
   const arma::umat   arma_perm = perm - 1;
   arma::uvec         row_index = arma::linspace<arma::uvec>(0, M - 1, M);
   arma::uvec         col_index(K);
   arma::vec          tmp(M);
   arma::vec          tmp2(K);
   arma::umat         index = arma::ones<arma::umat>(M, K);
   arma::umat         ind(M, K);
   arma::vec          dirich(K);
   arma::vec          shape(K);
   arma::vec          rate(K);
   arma::mat          func_val(M, K);

   for (unsigned int k = 0; k < K; ++k)
   {
      index.unsafe_col(k) *= k;
   }
   /* Set up the optimizer */
   std::vector<arma::mat*> f_data(2);

   f_data[0] = &lambda;
   f_data[1] = &weight;
   nlopt_opt optim;

   optim = nlopt_create(NLOPT_LN_NELDERMEAD, n);
   double lower_bound = 10e-6;

   nlopt_set_lower_bounds1(optim, lower_bound);
   nlopt_set_max_objective(optim, obj_stephens1997a_poisson, &f_data);
   while (value != value_next)
   {
      value = value_next;
      nlopt_optimize(optim, pars.memptr(), &value_next);
      for (unsigned int k = 0; k < K; ++k)
      {
         dirich.at(k) = pars[k];
         shape.at(k)  = pars[k + K];
         rate.at(k)   = pars[k + 2 * K];
      }
      /* Loop over permutations */
      for (unsigned int p = 0; p < P; ++p)
      {
         tmp                    = arma::prod(arma::exp(ldgamma(lambda(row_index, arma_perm.row(p)), shape, rate)), 1);
         tmp                   %= arma::exp(lddirichlet(weight(row_index, arma_perm.row(p)), dirich));
         func_val.unsafe_col(p) = arma::log(tmp);
      }
      for (unsigned int m = 0; m < M; ++m)
      {
         tmp2       = arma::conv_to<arma::vec>::from(func_val.row(m));
         col_index  = arma::sort_index(tmp2, "descend");
         ind.row(m) = arma_perm.row(col_index(0));
      }
      swapmat_by_index(lambda, ind);
      swapmat_by_index(weight, ind);
      swapumat_by_index(index, ind);
   }
   nlopt_destroy(optim);
   index += 1;
   return arma::conv_to<arma::imat>::from(index);
}


/**
 * ------------------------------------------------------------
 * stephens1997a_poisson_cc
 * @brief   Defines Stephens (1997a) relabelling algorithm for
 *          Binomial models. The nlopt library is used for
 *          optimization (Nelder-Mead algorithm)
 * @par values1 sampled lambda parameters; M x K
 * @par values2 sampled weight parameters; M x K
 * @par pars    Beta and Dirichlet hyper parameters
 * @par perm    matrix with all possible permutations of labels;
 * @return  matrix indicating the optimal labeling of sampled
 *          parameters; M x K
 * @detail  See Stephens (1997a)
 * @see nlopt
 * @author  Lars Simon Zehnder
 * ------------------------------------------------------------
 **/

// [[Rcpp::export]]

arma::imat stephens1997a_binomial_cc(Rcpp::NumericMatrix& values1,
                                     Rcpp::NumericMatrix values2,
                                     arma::vec pars, const arma::umat perm)
{
   const unsigned int M          = values1.rows();
   const unsigned int K          = values2.cols();
   const unsigned int P          = perm.n_rows;
   const unsigned int n          = pars.n_elem;
   double             value      = 1.0;
   double             value_next = -10.0e-8;
   arma::mat          pp(values1.begin(), M, K, true, true);
   arma::mat          weight(values2.begin(), M, K, true, true);
   const arma::umat   arma_perm = perm - 1;
   arma::uvec         row_index(M);
   arma::uvec         col_index(K);
   arma::vec          tmp(M);
   arma::vec          tmp2(K);
   arma::umat         index = arma::ones<arma::umat>(M, K);
   arma::umat         ind(M, K);
   arma::vec          dirich(K);
   arma::vec          shape1(K);
   arma::vec          shape2(K);
   arma::mat          func_val(M, K);

   for (unsigned int k = 0; k < K; ++k)
   {
      index.unsafe_col(k) *= k;
   }
   for (unsigned int m = 0; m < M; ++m)
   {
      row_index.at(m) = m;
   }
   /* Set up the optimizer */
   nlopt_opt optim;

   optim = nlopt_create(NLOPT_LN_NELDERMEAD, n);
   std::vector<arma::mat*> f_data(2);

   f_data[0] = &pp;
   f_data[1] = &weight;
   double lower_bound[1] = { 1e-10 };
   double upper_bound[1] = { 1e+7 };

   nlopt_set_lower_bounds(optim, lower_bound);
   nlopt_set_upper_bounds(optim, upper_bound);
   nlopt_set_max_objective(optim, obj_stephens1997a_binomial, &f_data);

   while (value != value_next)
   {
      value = value_next;
      nlopt_optimize(optim, pars.memptr(), &value_next);
      for (unsigned int k = 0; k < K; ++k)
      {
         dirich.at(k) = pars[k];
         shape1.at(k) = pars[k + K];
         shape2.at(k) = pars[k + 2 * K];
      }
      /* Loop over permutations */
      for (unsigned int p = 0; p < P; ++p)
      {
         tmp                    = arma::prod(arma::exp(ldbeta(pp(row_index, arma_perm.row(p)), shape1, shape2)), 1);
         tmp                   %= arma::exp(lddirichlet(weight(row_index, arma_perm.row(p)), dirich));
         func_val.unsafe_col(p) = arma::log(tmp);
      }
      for (unsigned int m = 0; m < M; ++m)
      {
         tmp2       = arma::conv_to<arma::vec>::from(func_val.row(m));
         col_index  = arma::sort_index(tmp2, "descend");
         ind.row(m) = arma_perm.row(col_index(0));
      }

      swapmat_by_index(pp, ind);
      swapmat_by_index(weight, ind);
      swapumat_by_index(index, ind);
   }
   nlopt_destroy(optim);
   index += 1;
   return arma::conv_to<arma::imat>::from(index);
}

// [[Rcpp::export]]

arma::imat stephens1997b_poisson_cc(Rcpp::NumericVector values,
                                    Rcpp::NumericMatrix comp_par,
                                    Rcpp::NumericMatrix weight_par,
                                    signed int max_iter = 200)
{
   unsigned int            N          = values.size();
   unsigned int            M          = comp_par.rows();
   unsigned int            K          = comp_par.cols();
   double                  value      = 1.0;
   double                  value_next = 0.0;
   arma::vec               arma_values(values.begin(), N, false, true);
   arma::mat               lambda(comp_par.begin(), M, K, true, true);
   arma::mat               weight(weight_par.begin(), M, K, true, true);
   arma::umat              index(M, K);
   arma::umat              index_out(M, K);
   arma::umat              indM(K, K);
   arma::mat               pmat_hat(N, K);
   arma::mat               cost(K, K);
   arma::uvec              seq_vec(K);
   std::vector<arma::mat*> mat_vector(M);

   pmat_hat  = arma::zeros(N, K);
   index_out = arma::ones<arma::umat>(M, K);
   for (unsigned int k = 0; k < K; ++k)
   {
      seq_vec.at(k)            = k * K;
      index_out.unsafe_col(k) *= (k + 1);
   }
   for (unsigned int m = 0; m < M; ++m)
   {
      arma::mat* pmat_ptr = new arma::mat(N, K);
      /* Save a pointer to the STL vector */
      mat_vector[m] = pmat_ptr;
   }
   signed int iter = 0;

   while (value != value_next)
   {
      iter      += 1;
      value      = value_next;
      value_next = 0.0;
      /* For all sampled MCMC parameters a matrix
       * pmat (N x K) is computed with p_ij
       * indicating the probability for a value i
       * being from component j.
       * */
      for (unsigned int m = 0; m < M; ++m)
      {
         for (unsigned int n = 0; n < N; ++n)
         {
            mat_vector[m]->row(n) = weight.row(m)
                                    % dpoisson(arma_values.at(n), lambda.row(m));
            mat_vector[m]->row(n) /= arma::sum(mat_vector[m]->row(n));
         }
         pmat_hat += *(mat_vector[m]);
      }
      /* This computes the reference estimator P_hat*/
      pmat_hat /= M;
      /* Now for each sampled MCMC parameter it is
       * searched for the optimal label by computing
       * the Kullback-Leibler distance of each 'pmat'
       * column 'l' from column 'k' of the reference
       * estimator P_hat.
       * The cost matrix cost_mat contains then the
       * distance of column 'l' from column 'k'.
       * An optimal assignment method computes the minimal
       * 'cost' regarding the labeling.
       * If 'k' is therein assigned to 'l', than the
       * label 'k' is switched to 'l'.
       * */
      for (unsigned int m = 0; m < M; ++m)
      {
         for (unsigned int k = 0; k < K; ++k)
         {
            for (unsigned int l = 0; l < K; ++l)
            {
               cost(k, l) = kulback_leibler(mat_vector[m]->unsafe_col(l),
                                            pmat_hat.unsafe_col(k));
            }
         }
         value_next += arma::trace(cost);
         /* Assignment */
         indM = hungarian(cost);
         arma::uvec f = arma::find(indM.t() == 1);
         index.row(m) = arma::trans(arma::find(indM.t() == 1) - seq_vec);
      }
      /* Permute parameters */
      swapmat_by_index(lambda, index);
      swapmat_by_index(weight, index);
      swapumat_by_index(index_out, index);
      pmat_hat.fill(0.0);
   }
   for (unsigned int m = 0; m < M; ++m)
   {
      delete mat_vector[m];
   }
   return arma::conv_to<arma::imat>::from(index_out);
}

// [[Rcpp::export]]

arma::imat stephens1997b_binomial_cc(Rcpp::NumericVector values,
                                     Rcpp::NumericVector reps, Rcpp::NumericMatrix comp_par,
                                     Rcpp::NumericMatrix weight_par)
{
   unsigned int            N          = values.size();
   unsigned int            M          = comp_par.rows();
   unsigned int            K          = comp_par.cols();
   double                  value      = 1.0;
   double                  value_next = 0.0;
   arma::vec               arma_values(values.begin(), N, false, true);
   arma::vec               arma_reps(reps.begin(), N, false, true);
   arma::mat               p(comp_par.begin(), M, K, true, true);
   arma::mat               weight(weight_par.begin(), M, K, true, true);
   arma::umat              index(M, K);
   arma::umat              index_out(M, K);
   arma::umat              indM(K, K);
   arma::mat               pmat_hat(N, K);
   arma::mat               cost(K, K);
   arma::uvec              seq_vec(K);
   std::vector<arma::mat*> mat_vector(M);

   pmat_hat  = arma::zeros(N, K);
   index_out = arma::ones<arma::umat>(M, K);
   for (unsigned int k = 0; k < K; ++k)
   {
      seq_vec.at(k)            = k * K;
      index_out.unsafe_col(k) *= (k + 1);
   }
   for (unsigned int m = 0; m < M; ++m)
   {
      arma::mat* pmat_ptr = new arma::mat(N, K);
      /* Save a pointer to the STL vector */
      mat_vector[m] = pmat_ptr;
   }
   while (value != value_next)
   {
      value      = value_next;
      value_next = 0.0;
      /* For all sampled MCMC parameters a matrix
       * pmat (N x K) is computed with p_ij
       * indicating the probability for a value i
       * being from component j.
       * */
      for (unsigned int m = 0; m < M; ++m)
      {
         for (unsigned int n = 0; n < N; ++n)
         {
            mat_vector[m]->row(n) = weight.row(m)
                                    % dbinomial(arma_values.at(n), arma_reps.at(n), p.row(m));
            mat_vector[m]->row(n) /= arma::sum(mat_vector[m]->row(n));
         }
      }
      for (unsigned int m = 0; m < M; ++m)
      {
         pmat_hat += *(mat_vector[m]);
      }
      /* This computes the reference estimator P_hat*/
      pmat_hat /= M;
      /* Now for each sampled MCMC parameter it is
       * searched for the optimal label by computing
       * the Kullback-Leibler distance of each 'pmat'
       * column 'l' from column 'k' of the reference
       * estimator P_hat.
       * The cost matrix cost_mat contains then the
       * distance of column 'l' from column 'k'.
       * An optimal assignment method computes the minimal
       * 'cost' regarding the labeling.
       * If 'k' is therein assigned to 'l', than the
       * label 'k' is switched to 'l'.
       * */
      for (unsigned int m = 0; m < M; ++m)
      {
         for (unsigned int k = 0; k < K; ++k)
         {
            for (unsigned int l = 0; l < K; ++l)
            {
               cost(k, l) = kulback_leibler(mat_vector[m]->unsafe_col(l),
                                            pmat_hat.unsafe_col(k));
            }
         }
         value_next += arma::trace(cost);
         /* Assignment */
         indM = hungarian(cost);
         arma::uvec f = arma::find(indM.t() == 1);
         index.row(m) = arma::trans(arma::find(indM.t() == 1) - seq_vec);
      }
      /* Permute parameters */
      swapmat_by_index(p, index);
      swapmat_by_index(weight, index);
      swapumat_by_index(index_out, index);
      pmat_hat.fill(0.0);
   }
   for (unsigned int m = 0; m < M; ++m)
   {
      delete mat_vector[m];
   }
   return arma::conv_to<arma::imat>::from(index_out);
}

// [[Rcpp::export]]

arma::imat stephens1997b_exponential_cc(Rcpp::NumericVector values,
                                        Rcpp::NumericMatrix comp_par,
                                        Rcpp::NumericMatrix weight_par)
{
   unsigned int            N          = values.size();
   unsigned int            M          = comp_par.rows();
   unsigned int            K          = comp_par.cols();
   double                  value      = 1.0;
   double                  value_next = 0.0;
   arma::vec               arma_values(values.begin(), N, false, true);
   arma::mat               lambda(comp_par.begin(), M, K, true, true);
   arma::mat               weight(weight_par.begin(), M, K, true, true);
   arma::umat              index(M, K);
   arma::umat              index_out(M, K);
   arma::umat              indM(K, K);
   arma::mat               pmat_hat(N, K);
   arma::mat               cost(K, K);
   arma::uvec              seq_vec(K);
   std::vector<arma::mat*> mat_vector(M);

   pmat_hat  = arma::zeros(N, K);
   index_out = arma::ones<arma::umat>(M, K);
   for (unsigned int k = 0; k < K; ++k)
   {
      seq_vec.at(k)            = k * K;
      index_out.unsafe_col(k) *= (k + 1);
   }
   for (unsigned int m = 0; m < M; ++m)
   {
      arma::mat* pmat_ptr = new arma::mat(N, K);
      /* Save a pointer to the STL vector */
      mat_vector[m] = pmat_ptr;
   }
   while (value != value_next)
   {
      value      = value_next;
      value_next = 0.0;
      /* For all sampled MCMC parameters a matrix
       * pmat (N x K) is computed with p_ij
       * indicating the probability for a value i
       * being from component j.
       * */
      for (unsigned int m = 0; m < M; ++m)
      {
         for (unsigned int n = 0; n < N; ++n)
         {
            mat_vector[m]->row(n) = weight.row(m)
                                    % dexponential(arma_values.at(n), lambda.row(m));
            mat_vector[m]->row(n) /= arma::sum(mat_vector[m]->row(n));
         }
      }
      for (unsigned int m = 0; m < M; ++m)
      {
         pmat_hat += *(mat_vector[m]);
      }
      /* This computes the reference estimator P_hat*/
      pmat_hat /= M;
      /* Now for each sampled MCMC parameter it is
       * searched for the optimal label by computing
       * the Kullback-Leibler distance of each 'pmat'
       * column 'l' from column 'k' of the reference
       * estimator P_hat.
       * The cost matrix cost_mat contains then the
       * distance of column 'l' from column 'k'.
       * An optimal assignment method computes the minimal
       * 'cost' regarding the labeling.
       * If 'k' is therein assigned to 'l', than the
       * label 'k' is switched to 'l'.
       * */
      for (unsigned int m = 0; m < M; ++m)
      {
         for (unsigned int k = 0; k < K; ++k)
         {
            for (unsigned int l = 0; l < K; ++l)
            {
               arma::vec mycol = mat_vector[m]->unsafe_col(l);
               cost(k, l) = kulback_leibler(mat_vector[m]->unsafe_col(l),
                                            pmat_hat.unsafe_col(k));
            }
         }
         value_next += arma::trace(cost);
         /* Assignment */
         indM = hungarian(cost);
         arma::uvec f = arma::find(indM.t() == 1);
         index.row(m) = arma::trans(arma::find(indM.t() == 1) - seq_vec);
      }
      /* Permute parameters */
      swapmat_by_index(lambda, index);
      swapmat_by_index(weight, index);
      swapumat_by_index(index_out, index);
      pmat_hat.fill(0.0);
   }
   for (unsigned int m = 0; m < M; ++m)
   {
      delete mat_vector[m];
   }
   return arma::conv_to<arma::imat>::from(index_out);
}

