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

#ifndef __FINMIX_OPTIMIZE_H__
#define __FINMIX_OPTIMIZE_H__

#include <vector>
#include "distributions.h"

// ============================================================
// Objective functions
// ------------------------------------------------------------

/**
 * ------------------------------------------------------------
 * @brief   Defines the objective function for the relabeling
 *          algorithm in Stephens (1997a) as prescribed by the
 *          nlopt library and to be applied on a Poisson model.
 * @par x       parameter vector to be optimized over
 * @par grad    gradient; not used
 * @par f_data  Armadillo matrix pointer to the data
 * @return  function value
 * @detail  The nlopt library needs an objective function with
 *          that returns a double value and uses certain
 *          parameters. In detail, the f_data pointer points to
 *          a std::vector filled with Armadillo matrix pointers.
 *          Inside the objective function a static_cast is used
 *          to access the data.
 *          The function calculates the sum over all logarithms
 *          of the Gamma prior added by the logarithms of the
 *          Dirichlet prior for the weights.
 * @see nlopt::max_objective
 * @author  Lars Simon Zehnder
 * ------------------------------------------------------------
 **/
inline
double obj_stephens1997a_poisson(unsigned n, const double* x,
                                 double* grad, void *f_data)
{
   std::vector<arma::mat*> *arma_data = static_cast<std::vector<arma::mat*>* >(f_data);
   const unsigned int      M          = (*arma_data)[0]->n_rows;
   const unsigned int      K          = (*arma_data)[0]->n_cols;
   arma::vec               rvalues(M);
   arma::vec               dirich(&x[0], K);
   arma::vec               shape(&x[0] + K, K);
   arma::vec               rate(&x[0] + 2 * K, K);

   rvalues  = lddirichlet((*(*arma_data)[1]), dirich);
   rvalues += arma::sum(ldgamma((*(*arma_data)[0]),
                                shape, rate), 1);
   if (rvalues.has_inf())
   {
      rvalues.elem(arma::find(rvalues == arma::datum::inf)).fill(10.0e+6);
      rvalues.elem(arma::find(rvalues == -arma::datum::inf)).fill(-10.0e+6);
   }
   else if (rvalues.has_nan())
   {
      rvalues.elem(arma::find(rvalues == arma::datum::nan)).zeros();
   }

   return arma::as_scalar(arma::sum(rvalues));
}
/**
 * ------------------------------------------------------------
 * @brief   Defines the objective function for the relabeling
 *          algorithm in Stephens (1997a) as prescribed by the
 *          nlopt library and to be applied on a Binomial model.
 * @par x       parameter vector to be optimized over
 * @par grad    gradient; not used
 * @par f_data  Armadillo matrix pointer to the data
 * @return  function value
 * @detail  The nlopt library needs an objective function with
 *          that returns a double value and uses certain
 *          parameters. In detail, the f_data pointer points to
 *          a std::vector filled with Armadillo matrix pointers.
 *          Inside the objective function a static_cast is used
 *          to access the data.
 *          The function calculates the sum over all logarithms
 *          of the Beta prior added by the logarithms of the
 *          Dirichlet prior for the weights.
 * @see nlopt::max_objective
 * @author  Lars Simon Zehnder
 * ------------------------------------------------------------
 **/
inline
double obj_stephens1997a_binomial(unsigned n, const double* x,
                                  double* grad, void *f_data)
{
   std::vector<arma::mat*> *arma_data = static_cast<std::vector<arma::mat*>* >(f_data);
   const unsigned int      M          = (*arma_data)[0]->n_rows;
   const unsigned int      K          = (*arma_data)[0]->n_cols;
   arma::vec               rvalues(M);
   arma::vec               arma_x(*x);
   arma::vec               dirich(&x[0], K);
   arma::vec               shape1(&x[0] + K, K);
   arma::vec               shape2(&x[0] + 2 * K, K);

   rvalues  = lddirichlet((*(*arma_data)[1]), dirich);
   rvalues += arma::sum(ldbeta((*(*arma_data)[0]),
                               shape1, shape2), 1);
   return arma::as_scalar(arma::sum(rvalues));
}

#endif /* __FINMIX_OPTIMIZE_H__ */



