/******************************************************************************
*
* TODO: Project Title
*
* Copyright (C) 2012-2013 Lars Simon Zehnder. All Rights Reserved.
* Web: -
*
* Author: Lars Simon Zehnder <simon.zehnder@gmail.com>
*
* This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
* WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
*
******************************************************************************/

#ifndef __FINMIX_ALGORITHMS_H__
#define __FINMIX_ALGORITHMS_H__

#include <RcppArmadillo.h>

inline
double kulback_leibler(const arma::vec &values, const arma::vec &base)
{
   const unsigned int N = values.n_elem;
   //const unsigned int K = values.n_elem;
   double             rvalue;

   if (arma::any(values == 0.0) || arma::any(base == 0.0))
   {
      arma::vec values_smoothed = (values * N + 1.0) / (N + 1.0);
      arma::vec base_smoothed   = (base * N + 1.0) / (N + 1.0);
      rvalue = arma::sum(values_smoothed % arma::log(values_smoothed / base_smoothed));
   }
   else
   {
      rvalue = arma::sum(values % arma::log(values / base));
   }

   //Rcpp::Rcout << "kullback: " << rvalue << std::endl;
   return rvalue;
}

inline
void swapmat_by_index(arma::mat &values, arma::umat index)
{
   const unsigned int K = values.n_cols;
   const unsigned int M = values.n_rows;
   arma::uvec         row_index(1);
   arma::urowvec      swap_index(K);

   for (unsigned int i = 0; i < M; ++i)
   {
      row_index.at(0) = i;
      swap_index      = index.row(i);
      values.row(i)   = values.submat(row_index, swap_index);
   }
}

inline
void swapumat_by_index(arma::umat &values, const arma::umat &index)
{
   const unsigned int K = values.n_cols;
   const unsigned int M = values.n_rows;
   arma::uvec         row_index(1);
   arma::urowvec      swap_index(K);

   for (unsigned int i = 0; i < M; ++i)
   {
      row_index.at(0) = i;
      swap_index      = index.row(i);
      values.row(i)   = values.submat(row_index, swap_index);
   }
}

#endif /* __FINMIX_ALGORITHMS_H__ */



