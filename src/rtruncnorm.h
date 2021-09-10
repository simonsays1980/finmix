/******************************************************************************
*
* TODO: Project Title
*
* Copyright (C) 2003-2009 ascolab GmbH. All Rights Reserved.
* Web: http://www.ascolab.com
*
* Author: Gerhard Gappmeier <gerhard.gappmeier@ascolab.com>
*
* This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
* WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
*
******************************************************************************/

#ifndef __FINMIX_RTRUNCNORM_H__
#define __FINMIX_RTRUNCNORM_H__

#define _USE_MATH_DEFINES
#define M1_SQRT_2PI    std::sqrt(2 * M_PI)

#include <cmath>
#include <RcppArmadillo.h>

static const double t1 = 0.15;
static const double t2 = 2.18;
static const double t3 = 0.725;
static const double t4 = 0.45;

inline
static double ers_a_inf(const double& a)
{
   const double ainv = 1.0 / a;
   double       x, rho;

   do
   {
      x   = R::rexp(ainv) + a;   /* rexp works with 1/lambda */
      rho = std::exp(-0.5 * std::pow((x - a), 2.0));
   } while (R::runif(0.0, 1.0) > rho);
   return x;
}

/* Exponential rejection sampling (a,b) */
inline
static double ers_a_b(const double& a, const double& b)
{
   const double ainv = 1.0 / a;
   double       x, rho;

   do
   {
      x   = R::rexp(ainv) + a; /* rexp works with 1/lambda */
      rho = exp(-0.5 * std::pow((x - a), 2.0));
   } while (R::runif(0.0, 1.0) > rho || x > b);
   return x;
}

/* Normal rejection sampling (a,b) */
inline
static double nrs_a_b(const double& a, const double& b)
{
   double x = a - 1.0;

   while (x < a || x > b)
   {
      x = R::rnorm(0.0, 1.0);
   }
   return x;
}

/* Normal rejection sampling (a,inf) */
inline
static double nrs_a_inf(const double & a)
{
   double x = a - 1.0;

   while (x < a)
   {
      x = R::rnorm(0.0, 1.0);
   }
   return x;
}

/* Half-normal rejection sampling */
inline
static double hnrs_a_b(const double& a, const double& b)
{
   double x = a - 1.0;

   while (x < a || x > b)
   {
      x = R::rnorm(0.0, 1.0);
      x = std::fabs(x);
   }
   return x;
}

/* Uniform rejection sampling */
inline
static double urs_a_b(const double& a, const double& b)
{
   const double phi_a = R::dnorm(a, 0.0, 1.0, 0);
   double       x     = 0.0;

   /* Upper bound of normal density on [a,b] */
   const double ub = a < 0.0 && b > 0.0 ? M1_SQRT_2PI : phi_a;

   do
   {
      x = R::runif(a, b);
   } while (R::runif(0.0, 1.0) * ub > R::dnorm(x, 0.0, 1.0, 0));
   return x;
}

/* Previously this was referred to as type 1 sampling: */
inline
static double r_lefttruncnorm(const double& a, const double& mean,
                              const double& sd)
{
   const double alpha = (a - mean) / sd;

   if (alpha < t4)
   {
      return mean + sd * nrs_a_inf(alpha);
   }
   else
   {
      return mean + sd * ers_a_inf(alpha);
   }
}

inline
static double r_righttruncnorm(const double& b, const double& mean,
                               const double& sd)
{
   const double beta = (b - mean) / sd;

   /* Exploit symmetry */
   return mean + sd * r_lefttruncnorm(-beta, 0.0, 1.0);
}

inline
static double r_truncnorm(const double& a, const double& b,
                          const double& mean, const double& sd)
{
   const double alpha = (a - mean) / sd;
   const double beta  = (b - mean) / sd;
   const double phi_a = R::dnorm(alpha, 0.0, 1.0, 0);
   const double phi_b = R::dnorm(beta, 0.0, 1.0, 0);

   if (beta <= alpha)
   {
      return NA_REAL;
   }
   else if (alpha <= 0.0 && 0 <= beta)
   {
      if (phi_a <= t1 || phi_b <= t1)
      {
         return mean + sd * nrs_a_b(alpha, beta);
      }
      else
      {
         return mean + sd * urs_a_b(alpha, beta);
      }
   }
   else if (alpha > 0)
   {
      if (phi_a / phi_b <= t2)
      {
         return mean + sd * urs_a_b(alpha, beta);
      }
      else
      {
         if (alpha < t3)
         {
            return mean + sd * hnrs_a_b(alpha, beta);
         }
         else
         {
            return mean + sd * ers_a_b(alpha, beta);
         }
      }
   }
   else
   {
      if (phi_b / phi_a <= t2)
      {
         return mean + sd * urs_a_b(-beta, -alpha);
      }
      else
      {
         if (beta > -t3)
         {
            return mean + sd * hnrs_a_b(-beta, -alpha);
         }
         else
         {
            return mean + sd * ers_a_b(-beta, -alpha);
         }
      }
   }
}

inline
arma::vec do_rtruncnorm(const unsigned int& n, const double& a,
                        const double& b, const double& mean, const double& sd)
{
   arma::vec      output(n);
   Rcpp::RNGScope scope;

   for (unsigned int i = 0; i < n; ++i)
   {
      if (R_FINITE(a) && R_FINITE(b))
      {
         output(i) = r_truncnorm(a, b, mean, sd);
      }
      else if (R_NegInf == a && R_FINITE(b))
      {
         output(i) = r_righttruncnorm(b, mean, sd);
      }
      else if (R_FINITE(a) && R_PosInf == b)
      {
         output(i) = r_lefttruncnorm(a, mean, sd);
      }
      else if (R_NegInf == a && R_PosInf == b)
      {
         output(i) = R::rnorm(mean, sd);
      }
      else
      {
         output(i) = NA_REAL;
      }
   }
   return output;
}

inline
double do_dtruncnorm(const double& x, const double& a,
                     const double& b, const double& mean,
                     const double& sd)
{
   double output = 0.0;

   if (a <= x && x <= b)    /* in range */
   {
      const double c1 = R::pnorm(a, mean, sd, 1, 0);
      const double c2 = R::pnorm(b, mean, sd, 1, 0);
      const double c3 = sd * (c2 - c1);
      const double c4 = R::dnorm((x - mean) / sd, 0.0, 1.0, 0);
      output = c4 / c3;
   }
   return output;
}
#endif /* __FINMIX_RTRUNCNORM_H__ */



