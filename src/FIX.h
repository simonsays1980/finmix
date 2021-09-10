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

#ifndef FIX_H
#define FIX_H

#include <RcppArmadillo.h>
#include "FinmixData.h"
#include "FinmixModel.h"
#include "FinmixPrior.h"
#include "FinmixMCMC.h"

// ===================================================================
// FIX mixin layer
// -------------------------------------------------------------------
/*
 * @brief   Mixin layer to implement the collaboration between 'Node'
 *          and 'Output' objects in case of Gibbs sampling with fixed
 *          indicators.
 * @par PriorType   parameter for prior distribution
 * @par ParType     parameter for posterior distribution
 * @par LogType     parameter for the log-likelihood function
 * @par ParOutType  parameter for the detailed storing of the parameters
 * @detail  Any implemented mixin layer describes the whole collabo-
 *          ration between 'Node' and 'Output' object to perform a
 *          Gibbs sampling of posterior parameters. The FIX mixin
 *          defines the two inner mixins 'Node' and 'Output' with
 *          variables needed to perform all actions for Gibbs
 *          sampling with fixed indicators (or for mixtures with one
 *          component only). These are e.g. variables needed to
 *          configure the algorithm, to perform random permutation
 *          Gibbs sampling, etc.
 *          The template parameters PriorType, ParType and LogType
 *          determine the specific model for that a Gibbs sampling
 *          should be performed. In particular they must specifiy
 *          parameters and an 'update()' function that can be called
 *          from the inner mixin 'Node's 'update()' function. The
 *          ParOutType parameter determines the specific storage
 *          prcocess for the parameters of a chosen model and has to be
 *          provided. It must contain a 'store()' method that can be
 *          called from the inner mixin 'Output's 'store()' function.
 * @see IND, HIER, POST, ADAPTER, BASE
 * @author Lars SImon Zehnder
 * ------------------------------------------------------------------
 */
template <typename PriorType, typename ParType, typename LogType,
          typename ParOutType>
class FIX {
public:
/**
 * ---------------------------------------------------------
 * Node mixin
 * ---------------------------------------------------------
 *
 * @brief   Holds all variables and method to perform the
 *          steps of a Gibbs sampler.
 * @detail  This class defines the variables needed for
 *          configuration of the algorithm as well as random
 *          permutation Gibbs sampling. The workhorse of this
 *          mixin is the virtual method 'update()' that
 *          performs the update step and calls any 'update()'
 *          function of related classes.
 * @see IND, HIER, POST, ADAPTER, BASE
 * --------------------------------------------------------
 */
class Node {
public:
const unsigned int K;
const unsigned int N;
const unsigned int M;
const unsigned int BURNIN;
const unsigned int STORES;
const bool INDICFIX;
const bool STARTPAR;
const bool HIER;
const bool RANPERM;
const bool STOREPOST;
PriorType hyperPar;
ParType par;
LogType log;
const arma::mat y;
arma::ivec S;
const arma::mat expos;
const arma::vec T;
arma::urowvec compIndex;
arma::urowvec permIndex;
arma::urowvec compIndex2;

Node (const FinmixData&, const FinmixModel&,
      const FinmixPrior&, const FinmixMCMC&);
virtual void update();
};
/**
 * -------------------------------------------------------
 * Output mixin
 * -------------------------------------------------------
 *
 * @brief   Stores all sampled parameters and additional
 *          information in container pointers.
 * @detail  This class defines container pointers needed
 *          to store any information from sampling.
 *          The workhorse of this inner mixin is the
 *          'store()' method that performs the storing
 *          process thereby calling all 'store()' methods
 *          of related classes.
 * @see IND, BASE, HIER, POST, ADAPTER
 * ------------------------------------------------------
 */
class Output {
public:
const unsigned int M;
const bool RANPERM;
ParOutType par;
arma::vec* mixlik;
arma::vec* mixprior;

Output (Rcpp::S4&);
virtual void store(const unsigned int&, Node&);
};
Node node;
Output output;

FIX (const FinmixData&, const FinmixModel&, const FinmixPrior&,
     const FinmixMCMC&, Rcpp::S4&);
virtual ~FIX ()
{
}
virtual void update();
virtual void store(const unsigned int&);
};

// ============================================================
// Node mixin definitions
// ------------------------------------------------------------

/**
 * ------------------------------------------------------------
 * Node::Node
 * @brief   Constructor of inner mixin 'Node'
 * @see FinmixModel, FinmixPrior, FinmixMCMC, HIER::Node::Node,
 *      IND::Node::Node, POST::Node::Node
 * @author  Lars Simon Zehnder
 * ------------------------------------------------------------
 **/
template <typename PriorType, typename ParType, typename LogType,
          typename ParOutType>
FIX <PriorType, ParType, LogType, ParOutType>::Node::Node (const FinmixData& data,
                                                           const FinmixModel& model, const FinmixPrior& prior,
                                                           const FinmixMCMC& mcmc) :
   K(model.K), N(data.N), M(mcmc.M), BURNIN(mcmc.burnIn),
   STORES(mcmc.storeS), INDICFIX(model.indicFix),
   STARTPAR(mcmc.startPar), HIER(prior.hier),
   RANPERM(mcmc.ranPerm), STOREPOST(mcmc.storePost),
   hyperPar(prior), par(mcmc.startPar, model), log(),
   y(data.y), S(data.S), expos(data.expos), T(data.T),
   compIndex(model.K), permIndex(model.K), compIndex2(model.K)
{
   for (unsigned int k = 0; k < K; ++k)
   {
      compIndex(k) = k;
   }
}

/**
 * -----------------------------------------------------------
 * Node::update
 * @brief   Updates the 'node' object.
 * @detail  Virtual. Performs any updates on parameters and
 *          then starts random permutation of sampled parameters.
 *          It is this function, that is passed forward via
 *          inheritance to any other mixin.
 * @see IND::Node::update, POST::Node::update,
 *      HIER::Node::update
 * @author  Lars Simon Zehnder
 * -----------------------------------------------------------
 **/
template <typename PriorType, typename ParType, typename LogType,
          typename ParOutType>
void FIX <PriorType, ParType, LogType, ParOutType>::Node::update()
{
   hyperPar.update(K, y, S, T, par);
   par.update(hyperPar);
   hyperPar.updateHier(par);
   log.update(K, y, S, expos, T, par, hyperPar);
   if (RANPERM && K > 1)
   {
      permIndex  = arma::shuffle(compIndex, 1);
      compIndex2 = (permIndex == compIndex);
      if (arma::sum(compIndex) != K)
      {
         par.permute(compIndex, permIndex);
      }
   }
}

// ==========================================================
// Output mixin definitions
// ----------------------------------------------------------

/**
 * ----------------------------------------------------------
 * Output::Output
 * @brief   Constructs an object of class 'Output' inside of
 *          the mixin layer.
 * @par classS4     object of class Rcpp::S4
 * @detail  'classS4' is an R S4 class object wrapped by an
 *          Rcpp::S4 object holding a certain structure of
 *          containers to store sampled parameters, log-like-
 *          lihoods, etc. Note, the Rcpp::S4 object references
 *          in its objects to memory allocated in R. To avoid
 *          copying memory, pointers are used to represent the
 *          containers in the C++ application. For each
 *          Armadillo object its advanced constructor is
 *          called to reuse auxiliary memory and fix the size.
 * @see FIX::Output::Output, HIER::Output::Output,
 *      POST::Output::Output, Rcpp::S4, ?S4 (in R), arma::mat
 * @author  Lars Simon Zehnder
 * ----------------------------------------------------------
 **/
template <typename PriorType, typename ParType, typename LogType, typename ParOutType>
FIX <PriorType, ParType, LogType, ParOutType>::Output::Output (Rcpp::S4& classS4) :
   M(Rcpp::as<unsigned int>((SEXP)classS4.slot("M"))),
   RANPERM(Rcpp::as<bool>((SEXP)classS4.slot("ranperm"))),
   par(Rcpp::as<Rcpp::List>((SEXP)classS4.slot("par")))
{
   Rcpp::List          tmpLog((SEXP)classS4.slot("log"));
   Rcpp::NumericVector tmpMixLik((SEXP)tmpLog["mixlik"]);
   Rcpp::NumericVector tmpMixPrior((SEXP)tmpLog["mixprior"]);

   mixlik   = new arma::vec(tmpMixLik.begin(), M, false, true);
   mixprior = new arma::vec(tmpMixPrior.begin(), M, false, true);
}

/**
 * ---------------------------------------------------------
 * Output::store
 * @brief   Stores the sampled parameters into containers.
 * @par m       iteration count
 * @par node    object this->Node
 * @detail  Takes the iteration number and a 'Node' object
 *          holding all information from one sampling step
 *          and stores it to the containers pointed to in-
 *          side the 'Output' class. It thereby always
 *          checks if the iteration is part of the burnin
 *          phase or the sampling phase, if indicators
 *          should be stored at all, etc.
 * @see IND::Output::store, HIER::Output::store,
 *      POST::Output::store,
 * @author Lars Simon Zehnder
 * ---------------------------------------------------------
 **/
template <typename PriorType, typename ParType, typename LogType,
          typename ParOutType>
void FIX <PriorType, ParType, LogType, ParOutType>::Output::store(const unsigned
                                                                  int& m, Node& node)
{
   if (m >= node.BURNIN)
   {
      const unsigned int index = m - node.BURNIN;
      (*mixlik)(index)   = node.log.mixlik;
      (*mixprior)(index) = node.log.mixprior;
      par.store(index, node.par);
   }
}

// ========================================================
// FIX mixin layer
// --------------------------------------------------------

/**
 * --------------------------------------------------------
 * FIX<Super>::FIX
 * @brief   Constructs an object of the parameterized mixin
 *          layer.
 * @par data    object of class FinmixData, holds the data
 * @par model   object of class FinmixModel, holds model
 *              information
 * @par prior   object of class FinmixPrior, holds prior
 *              information
 * @par mcmc    object of class FinmixMCMC, holds info for
 *              algorithmic configurations
 * @par classS4 object of class Rcpp::S4 to pass output
 *              container pointer
 * @detail  Note, that this constructor must include all
 *          parameters needed in construction of the inner
 *          mixins.
 * @see FinmixData, FinmixModel, FinmixPrior,
 *      FinmixMCMC, Rcpp::S4, IND<Super>::IND,
 *      POST<Super>::POST, HIER<Super>::HIER
 * @author  Lars Simon Zehnder
 * -------------------------------------------------------
 **/
template <typename PriorType, typename ParType, typename LogType,
          typename ParOutType>
FIX <PriorType, ParType, LogType, ParOutType>::FIX (const FinmixData& data,
                                                    const FinmixModel& model, const FinmixPrior& prior, const FinmixMCMC&
                                                    mcmc, Rcpp::S4& classS4) :
   node(data, model, prior, mcmc), output(classS4)
{
}

/**
 * -------------------------------------------------------
 * IND<Super>::update
 * @brief   Triggers the update process for each step
 *          the sampler. Passes responsibility to 'Node's
 *          'update()' method.
 * @see Node::update
 * @author  Lars Simon Zehnder
 * -------------------------------------------------------
 **/
template <typename PriorType, typename ParType, typename LogType,
          typename ParOutType>
void FIX <PriorType, ParType, LogType, ParOutType>::update()
{
   node.update();
}

/**
 * -------------------------------------------------------
 * IND<Super>::store
 * @brief   Triggers the store process for each step of
 *          the sampler. Passes responsibility to 'Output's
 *          'store()' method.
 * @see Output::store
 * @author  Lars Simon Zehnder
 * -------------------------------------------------------
 **/
template <typename PriorType, typename ParType, typename LogType,
          typename ParOutType>
void FIX <PriorType, ParType, LogType, ParOutType>::store(const unsigned int& m)
{
   output.store(m, node);
}
#endif
