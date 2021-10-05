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
#include "distributions.h"
#include "hungarian.h"
#include "mincol.h"
#include "moments.h"

//' Swaps values in each row
//' 
//' @description
//' This function swaps the values in each row of a matrix by permuting the 
//' columns via the indices provided in the `index` matrix. All 
//' `swapElements()`-methods use this function internally. The code is extended 
//' to `C++` using the packages `Rcpp` and `RcppArmadillo`.
//' 
//' @param values A matrix containing the values to be swapped. 
//' @param index An integer matrix defining how values should be swapped. 
//' @return A matrix with swapped values. 
//' @export
//' 
//' @examples 
//' values <- matrix(rnorm(10), nrow = 2)
//' index <- matrix(c(2,1), nrow = 5, ncol = 2)
//' swap_cc(values, index)
//' 
//' @seealso
//' * [swapElements()][mcmcoutput_class] for the calling function
// [[Rcpp::export]]
Rcpp::NumericMatrix swap_cc(Rcpp::NumericMatrix values, Rcpp::IntegerMatrix index)
{
   /* If dimensions of both arguments do not agree throw an exception */
   if (values.nrow() != index.nrow() || values.ncol() != index.ncol())
   {
      throw Rcpp::exception("Matrix dimensions disagree.");
   }
   /* Do not reuse memory from R as otherwise existing objects
    * get manipulated */
   const unsigned int K = values.ncol();
   const unsigned int M = values.nrow();
   arma::mat          values_arma(values.begin(), M, K, true, true);
   arma::imat         index_arma(index.begin(), M, K, true, true);
   arma::mat          values_copy(M, K);
   arma::umat         index_umat = arma::conv_to<arma::umat>::from(index_arma) - 1;
   arma::uvec         row_index(1);
   arma::urowvec      swap_index(K);

   for (unsigned int i = 0; i < M; ++i)
   {
      row_index.at(0)    = i;
      swap_index         = index_umat.row(i);
      values_copy.row(i) =
         values_arma.submat(row_index, swap_index);
   }
   return Rcpp::wrap(values_copy);
}

//' Swap elements in a 3d array
//' 
//' @description
//' This function swaps the elements in a three-dimensional array by using the 
//' scheme provided in the `index` matrix. 
//' 
//' @param values An array of dimension `M x r x K` of values to swap. 
//' @param index An integer matrix of dimension `M x K`. containing the scheme 
//'   by which values should be swapped.
//' @param A three-dimensional array with swapped values.
//' @export
//' 
//' @examples
//' values <- array(rnorm(40), dim = c(10, 2, 2))
//' index <- matrix(c(1,2), nrow = 10, ncol = 2)
//' swap_3d_cc(values, index)
//' 
//' @seealso 
//' * [swapElements()][mcmcoutput_class] for the calling method
//' * [swap_cc()] for the equivalent function for 2-dimensional arrays 
// [[Rcpp::export]]
Rcpp::NumericVector swap_3d_cc(Rcpp::NumericVector values, Rcpp::IntegerMatrix index)
{
   Rcpp::IntegerVector valDim = values.attr("dim");
   const unsigned int  M      = valDim[0];
   const unsigned int  r      = valDim[1];
   const unsigned int  K      = valDim[2];

   /* If dimensions of both arguments do not agree throw an exception */
   if (M != (unsigned)index.nrow() || K != (unsigned)index.ncol())
   {
      throw Rcpp::exception("Matrix dimensions disagree.");
   }
   arma::cube values_arma(values.begin(), M, r, K, false, true);
   arma::imat index_arma(index.begin(), M, K, false, true);
   arma::cube output(M, r, K);

   output.fill(0.0);
   arma::umat  index_umat = arma::conv_to<arma::umat>::from(index_arma) - 1;
   arma::umat  ik(M, K);
   arma::ucube ikr(M, 1, K);
   arma::ucube ikr2(M, r, K);
   arma::cube  ikr3(M, r, K);

   for (unsigned int k = 0; k < K; ++k)
   {
      ik                   = (index_arma - 1) == k;
      ikr.slices(0, K - 1) = ik;
      ikr2                 = arma::resize(ikr, M, r, K);
      for (unsigned int rr = 1; rr < r; ++rr)
      {
         ikr2.tube(0, rr, M - 1, rr) = ikr2.tube(0, 0, M - 1, 0);
      }
      ikr3  = arma::conv_to<arma::cube>::from(ikr2);
      ikr3 %= values_arma;
      for (unsigned int l = 0; l < K; ++l)
      {
         output.slice(k) += ikr3.slice(l);
      }
      ik.fill(0);
      ikr.fill(0);
      ikr2.fill(0);
      ikr3.fill(0);
   }
   return Rcpp::wrap(output);
}

//' Swap values in an integer matrix
//' 
//' @description
//' This function swaps the values in an integer matrix column-wise defined 
//' by the `index` matrix. This function is used mainly for the 
//' `swapElements()`-method of MCMC samples to swap the indicator values.
//' 
//' @param values An integer matrix containing the values to swap. 
//' @param index An integer matrix containing the indices by which values 
//'   should be swapped.
//' @return An integer matrix containing the swapped values.
//' @export
//' 
//' @examples 
//' values <- matrix(c(2, 4, 1, 3), nrow = 10, ncol = 2)
//' index <- matrix(c(1, 2), nrow = 10, ncol = 2)
//' swapInteger_cc(values, index)
//' 
//' @seealso 
//' * [swap_cc()] for the equivalent function for numeric values
//' * [swap_3d_cc()] for the equivalent function for three-dimensional arrays
// [[Rcpp::export]]
Rcpp::IntegerMatrix swapInteger_cc(Rcpp::IntegerMatrix values, Rcpp::IntegerMatrix index)
{
   /* If dimensions of both arguments do not agree throw an exception */
   if (values.nrow() != index.nrow() || values.ncol() != index.ncol())
   {
      throw Rcpp::exception("Matrix dimensions disagree.");
   }
   /* Do not reuse memory from R as otherwise existing objects
    * get manipulated */
   const unsigned int K = values.ncol();
   const unsigned int M = values.nrow();
   arma::imat         values_arma(values.begin(), M, K, true, true);
   arma::imat         index_arma(index.begin(), M, K, true, true);
   arma::imat         values_copy(M, K);
   arma::umat         index_umat = arma::conv_to<arma::umat>::from(index_arma) - 1;
   arma::uvec         row_index(1);
   arma::urowvec      swap_index(K);

   for (unsigned int i = 0; i < M; ++i)
   {
      row_index.at(0)    = i;
      swap_index         = index_umat.row(i);
      values_copy.row(i) =
         values_arma.submat(row_index, swap_index);
   }
   return Rcpp::wrap(values_copy);
}

//' Swap values of stored indicators
//' 
//' @description
//' This function is used to swap elements in the stored indicators from MCMC 
//' sampling. Note that this function reuses R memory and should therefore be 
//' treated with caution. Do not use this function unless you really know what 
//' you are doing. 
//' 
//' @param values An integer matrix containing the last indicators stored in 
//'   MCMC sampling. The number of these last stored indicators is defined by 
//'   the hpyer-parameter `storeS` in the `mcmc` object.
//' @param index An integer matrix defining the swapping scheme. 
//' @return A matrix with swapped values.
//' @export
//' 
//' @seealso
//' * [mcmc()] for the hyper-parameter `storeS`
//' * [swapElements()][mcmcoutput_class] for the calling method 
//' * [swapInteger_cc()] for the equivalent function that swaps simple integer 
//'   matrices
//' * [swap_3d_cc()] for a function that swaps values in three-dimensional 
//'   arrays
// [[Rcpp::export]]
Rcpp::IntegerMatrix swapInd_cc(Rcpp::IntegerMatrix values, Rcpp::IntegerMatrix index)
{
   /* If dimensions of both arguments do not agree throw an exception */
   if (values.ncol() != index.nrow())
   {
      throw Rcpp::exception("Matrix dimensions disagree.");
   }
   /* Reuse memory from R */
   const unsigned int N      = values.nrow();
   const unsigned int STORES = values.ncol();
   const unsigned int M      = index.nrow();
   const unsigned int K      = index.ncol();
   arma::imat         values_arma(values.begin(), N, STORES, true, true);
   arma::imat         index_arma(index.begin(), M, K, true, true);
   arma::imat         values_copy(N, STORES);

   for (unsigned int s = 0; s < STORES; ++s)
   {
      for (unsigned int i = 0; i < N; ++i)
      {
         values_copy(i, s) = (int)index_arma(s, (unsigned int)
                                             values_arma(i, s) - 1);
      }
   }
   return Rcpp::wrap(values_copy);
}

//' Swap the `ST` slot in the MCMC output
//' 
//' @description
//' This function is used to swap the elements in slot `ST` of an `mcmcoutput` 
//' object (An MCMC sampling output). The main difference to the 
//' [swapInteger_cc()] function is that this function reuses memory from R. Do 
//' only use this function, if you really know what you are doing.
//' 
//' @param values An integer matrix containing the values to swap in R memory.
//' @param index An integer matrix containing the swapping scheme. 
//' @return An integer matrix with swapped values.
//' @export
//' 
//' @seealso 
//' * [swapInteger_cc()] for the equivalent function not using R memory
//' * [swap_3d_cc()] for an equivalent function for three-dimensional arrays 
//' * [swapElements()][mcmcoutput_class] for the calling method
// [[Rcpp::export]]
Rcpp::IntegerVector swapST_cc(Rcpp::IntegerVector values, Rcpp::IntegerMatrix index)
{
   /* If dimensions of both arguments do not agree throw an exception */
   if (values.size() != index.nrow())
   {
      throw Rcpp::exception("Matrix dimensions disagree.");
   }
   /* Reuse memory from R */
   const unsigned int M = values.size();
   const unsigned int K = index.ncol();
   arma::ivec         values_arma(values.begin(), M, false, true);
   arma::imat         index_arma(index.begin(), M, K, false, true);
   arma::ivec         values_copy(M);

   for (unsigned int i = 0; i < M; ++i)
   {
      values_copy(i) = index_arma(i, (unsigned int)
                                  values_arma(i) - 1);
   }
   return Rcpp::wrap(values_copy);
}

//' Computes the log density of the Gamma distribution 
//' 
//' @description
//' For each shape and rate parameter pair the log gamma density is computed. 
//' Inside the function the unsafe access functions of Armadillo `at()` and 
//' `unsafe_col()` are used, so now boundary check is performed. In each step 
//' the `lngamma()` function from Rcpp's `R` namespace is used. At this time 
//' unused.
//' 
//' @param values A matrix of dimension `M x K` for which the log-density 
//'   should be calculated. 
//' @param shape A vector of dimension `K x 1` with Gamma shape parameters.
//' @param rate A vector of dimension `K x 1` with Gamma rate parameters.
//' @return A matrix of Gamma log-density values for each pair of parameters 
//'   in a column.
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix ldgamma_cc(Rcpp::NumericMatrix values, Rcpp::NumericVector shape,
                               Rcpp::NumericVector rate)
{
   /* Reuse memory from R */
   const unsigned int M = values.nrow();
   const unsigned int K = values.ncol();
   arma::mat          arma_values(values.begin(), M, K, false, true);
   arma::vec          arma_shape(shape.begin(), K, false, true);
   arma::vec          arma_rate(rate.begin(), K, false, true);
   arma::mat          arma_return(M, K);

   arma_return = ldgamma(arma_values, arma_shape, arma_rate);
   return Rcpp::wrap(arma_return);
}

//' Computes the density of the Gamma distribution 
//' 
//' @description
//' For each shape and rate parameter pair the gamma density is computed. 
//' Inside the function the unsafe access functions of Armadillo `at()` and 
//' `unsafe_col()` are used, so now boundary check is performed. In each step 
//' the `lngamma()` function from Rcpp's `R` namespace is used. At this time 
//' unused.
//' 
//' @param values A matrix of dimension `M x K` for which the density 
//'   should be calculated. 
//' @param shape A vector of dimension `K x 1` with Gamma shape parameters.
//' @param rate A vector of dimension `K x 1` with Gamma rate parameters.
//' @return A matrix of Gamma density values for each pair of parameters 
//'   in a column.
//' @export
// [[Rcpp::export]]
arma::mat dgamma_cc(Rcpp::NumericMatrix values, Rcpp::NumericVector shape,
                    Rcpp::NumericVector rate)
{
   /* Reuse memory from R */
   const unsigned int M = values.nrow();
   const unsigned int K = values.ncol();
   arma::mat          arma_values(values.begin(), M, K, false, true);
   arma::vec          arma_shape(shape.begin(), K, false, true);
   arma::vec          arma_rate(rate.begin(), K, false, true);
   arma::mat          arma_return(M, K);

   arma_return = exp(ldgamma(arma_values, arma_shape, arma_rate));
   return arma_return;
}

//' Computes the log density of the Dirichlet distribution 
//' 
//' @description
//' For each shape and rate parameter pair the log-Dirichlet density is 
//' computed. Inside the function the unsafe access functions of Armadillo 
//' `at()` and `unsafe_col()` are used, so now boundary check is performed. 
//' In each step the `lgammafn()` function from Rcpp's `R` namespace is used. 
//' At this time unused.
//' 
//' @param values A matrix of dimension `M x K` for which the log-density 
//'   should be calculated. 
//' @param par A vector of dimension `K x 1` containing the Dirichlet 
//'   parameters.
//' @return A vector of Dirichlet log-density values. 
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector lddirichlet_cc(Rcpp::NumericMatrix values, Rcpp::NumericVector par)
{
   /* Reuse memory from R */
   const unsigned int M = values.nrow();
   const unsigned int K = values.ncol();
   arma::mat          arma_values(values.begin(), M, K, false, true);
   arma::vec          arma_par(par.begin(), K, false, true);
   arma::vec          arma_return(M);

   arma_return = lddirichlet(arma_values, arma_par);
   return Rcpp::wrap(arma_return);
}

//' Computes the density of the Dirichlet distribution 
//' 
//' @description
//' For each shape and rate parameter pair the Dirichlet density is 
//' computed. Inside the function the unsafe access functions of Armadillo 
//' `at()` and `unsafe_col()` are used, so now boundary check is performed. 
//' In each step the `lgammafn()` function from Rcpp's `R` namespace is used. 
//' At this time unused.
//' 
//' @param values A matrix of dimension `M x K` for which the log-density 
//'   should be calculated. 
//' @param par A vector of dimension `K x 1` containing the Dirichlet 
//'   parameters.
//' @return A vector of Dirichlet density values. 
//' @export
// [[Rcpp::export]]
arma::vec ddirichlet_cc(Rcpp::NumericMatrix values, Rcpp::NumericVector par)
{
   /* Reuse memory from R */
   const unsigned int M = values.nrow();
   const unsigned int K = values.ncol();
   arma::mat          arma_values(values.begin(), M, K, false, true);
   arma::vec          arma_par(par.begin(), K, false, true);
   arma::vec          arma_return(M);

   arma_return = arma::exp(lddirichlet(arma_values, arma_par));
   return arma_return;
}

//' Compute the hungarian matrix
//' 
//' @description
//' This function calls an implementation of the Hungarian algorithm by Munkres. 
//' The Hungarian algorithm solves a weighted assignment problem on a bipartite 
//' graph. Note, here this algorithm is used in the re-labeling algorithm by 
//' Stephens (1997b).
//' 
//' @param cost A matrix containing the costs for each row source and column 
//'   target. 
//' @return An integer matrix defining the best solution to the assignment 
//'   problem.
//' @export
//' @seealso 
//' * [mcmcpermute()] for the calling function 
//' * [mcmcestimate()] for the function that uses the re-labeling algorithm by 
//'   Stephens (1997b)
//' 
//' @references
//' * Stephens, Matthew (1997b), "Dealing with Label-Switching in Mixture 
//'   Models", Journal of the Royal Statistical Society Series B, 62(4)
// [[Rcpp::export]]
arma::imat hungarian_cc(const arma::mat cost)
{
   arma::umat indM = hungarian(cost);

   return arma::conv_to<arma::imat>::from(indM);
}

//' Calculate moments on samples of multivariate mixture models 
//' 
//' @description
//' This function calculates the moments for MCMC samples of multivariate 
//' mixture models. Moments like means, standard deviations, kurtosis and 
//' skewness are computed for each iteration in MCMC sampling. The moments are 
//' used when plotting the traces of an MCMC sample output. 
//' 
//' @param classS4 An `mcmcoutput` class containing the MCMC samples.
//' @return A named list with vectors containing the data moments for each 
//'   iteration in the MCMC sample.
//' @export
//' @seealso 
//' * [mcmcoutput][mcmcoutput_class] for the `mcmcoutput` class definition
//' * [mixturemcmc()] for performing MCMC sampling
//' * [plotTraces][mcmcoutput_class] for the calling function
// [[Rcpp::export]]
Rcpp::List moments_cc(Rcpp::S4 classS4)
{
   Rcpp::S4   model    = Rcpp::as<Rcpp::S4>((SEXP)classS4.slot("model"));
   const bool indicfix = Rcpp::as<bool>((SEXP)model.slot("indicfix"));

   if (indicfix)
   {
      return Rcpp::wrap(moments_fix_cc(classS4));
   }
   else
   {
      return Rcpp::wrap(moments_ind_cc(classS4));
   }
}

//' Calculate moments on permuted samples of multivariate mixture models 
//' 
//' @description
//' This function calculates the moments for re-labeled MCMC samples of 
//' multivariate mixture models. Moments like means, standard deviations, 
//' kurtosis and skewness are computed for each iteration in MCMC sampling. The 
//' moments are used when plotting the traces of an MCMC sample output. 
//' 
//' @param classS4 An `mcmcoutputperm` class containing the re-labeled MCMC 
//'   samples.
//' @return A named list with vectors containing the data moments for each 
//'   iteration in the re-labeled MCMC sample.
//' @export
//' @seealso 
//' * [mcmcoutputperm][mcmcoutputperm_class] for the `mcmcoutput` class definition
//' * [mixturemcmc()] for performing MCMC sampling
//' * [mcmcpermute()] for re-labeling MCMC samples
//' * [plotTraces][mcmcoutputperm_class] for the calling function
// [[Rcpp::export]]
Rcpp::List permmoments_cc(Rcpp::S4 classS4)
{
   Rcpp::S4   model    = Rcpp::as<Rcpp::S4>((SEXP)classS4.slot("model"));
   const bool indicfix = Rcpp::as<bool>((SEXP)model.slot("indicfix"));

   if (indicfix)
   {
      return Rcpp::wrap(permmoments_fix_cc(classS4));
   }
   else
   {
      return Rcpp::wrap(permmoments_ind_cc(classS4));
   }
}