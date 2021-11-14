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
#ifndef POST_H
#define POST_H

#include <RcppArmadillo.h>

// ===================================================================
// POST mixin layer
// -------------------------------------------------------------------
/*
 * @brief   Mixin layer to implement the collaboration between 'Node'
 *          and 'Output' objects in case of Gibbs sampling when
 *          posterior hyper parameters should be stored.
 * @par Super   next mixin layer in the application
 * @detail  Any implemented mixin layer describes the whole collabo-
 *          ration between 'Node' and 'Output' object to perform a
 *          Gibbs sampling with storage of posterior hyper parameters.
 *          The mixin layer refines thereby its 'Super' class by
 *          defining new inner mixins 'Node' and 'Output' with
 *          additional variables needed to perform all actions for
 *          Gibbs sampling with storage of posterior hyper parameters.
 *          These are e.g. the additional container in the 'Output'
 *          mixin to store posterior hyper parameters.
 * @see FIX, IND, HIER, ADAPTER, BASE
 * @author Lars SImon Zehnder
 * ------------------------------------------------------------------
 */
template <typename Super, typename PostOutType>
class POST : public Super {
public:
/**
 * ---------------------------------------------------------
 * Node mixin
 * ---------------------------------------------------------
 *
 * @brief   Holds all variables and method to perform the
 *          steps of a Gibbs sampler.
 * @detail  This class inherits directly from the 'Super'
 *          class' 'Node' mixin. The workhorse of this mixin
 *          is the inherited virtual method 'update()' that
 *          performs the update step and calls any 'update()'
 *          function of related classes. Hierarchical
 *          parameters are then updated in the related
 *          classes whose 'update()' method gets called and
 *          knows what to do.
 * @see FIX, IND, POST, ADAPTER, BASE
 * --------------------------------------------------------
 */
class Node : public Super::Node {
public:
Node (const FinmixData&,
      const FinmixModel&,
      const FinmixPrior&,
      const FinmixMCMC&);
virtual ~Node ()
{
}
};
/**
 * -------------------------------------------------------
 * Output mixin
 * -------------------------------------------------------
 *
 * @brief   Stores all sampled parameters and additional
 *          information in container pointers.
 * @detail  This class inherits directly from the 'Super'
 *          class' 'Output' mixin. It defines the new
 *          container pointers needed to store any addi-
 *          tional information for sampling hyper
 *          parameters.Reusable functionality is inherited
 *          from 'Super's 'Output' class. The workhorse of
 *          this inner mixin is the 'store()' method that
 *          performs the storing process thereby calling
 *          all 'store()' methods of related classes.
 * @see FIX, BASE, IND, HIER, ADAPTER
 * ------------------------------------------------------
 */
class Output : public Super::Output {
public:
PostOutType post;

Output (Rcpp::S4&);
virtual ~Output ()
{
}
virtual void store(const
                   unsigned int&,
                   Node&);
};
Node node;
Output output;

POST (const FinmixData&, const FinmixModel&,
      const FinmixPrior&, const FinmixMCMC&,
      Rcpp::S4&);
virtual ~POST ()
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
 * @detail  Calls in its initialization list the constructor of
 *          its super class that takes the same parameters.
 * @see Super::Node::Node, FinmixModel, FinmixPrior,
 *      FinmixMCMC
 * @author  Lars Simon Zehnder
 * ------------------------------------------------------------
 **/
template <typename Super, typename PostOutType>
POST <Super, PostOutType>::Node::Node (const FinmixData& data,
                                       const FinmixModel& model, const FinmixPrior& prior,
                                       const FinmixMCMC& mcmc) :
   Super::Node(data, model, prior, mcmc)
{
}

/**
 * -----------------------------------------------------------
 * Node::update
 * @brief   Updates the 'node' object.
 * @detail  Virtual. Calls 'Super::Node::update()
 * @see Super::Node::update()
 * @author  Lars Simon Zehnder
 * -----------------------------------------------------------
 **/

// ==========================================================
// Output mixin definitions
// ----------------------------------------------------------

/**
 * ----------------------------------------------------------
 * Output::Output
 * @brief   Constructs an object of class 'Output' inside of
 *          the mixin layer.
 * @par classS4 object of class Rcpp::S4
 * @detail  Calls in its initialization list the constructor
 *          of the super class that takes the same parameter.
 *          'classS4' is an R S4 class object holding a
 *          certain structure of containers to store sampled
 *          parameters, log-likelihoods, etc. Note, the
 *          Rcpp::S4 object references in its objects to
 *          memory allocated in R. To avoid copying pointers
 *          are used to represent the containers in the C++
 *          application. For each Armadillo object its
 *          advanced constructor is called to reuse auxiliary
 *          memory and fix the size. For hierarchical prior
 *          modeling an additional container 'hyper' has to
 *          be prepared to store posterior hyper parameters.
 * @see FIX::Output::Output, IND::Output::Output,
 *      HIER::Output::Output, Rcpp::S4, ?S4 (in R),
 *      arma::mat
 * @author  Lars Simon Zehnder
 * ----------------------------------------------------------
 **/
template <typename Super, typename PostOutType>
POST <Super, PostOutType>::Output::Output (Rcpp::S4& classS4) :
   Super::Output(classS4),
   post(Rcpp::as<Rcpp::List>((SEXP)classS4.slot("post")))
{
}

/**
 * ---------------------------------------------------------
 * Output::store
 * @brief   Stores the sampled parameters into containers.
 * @par m       iteration count
 * @par node    object of class this->Node
 * @detail  Takes the iteration number and a 'Node' object
 *          holding all information from one sampling step
 *          and stores it to the containers pointed to in-
 *          side the 'Output' class. It thereby always
 *          checks if the iteration is part of the burnin
 *          phase or the sampling phase, etc.
 * @see FIX::Output::store, IND::Output::store,
 *      HIER::Output::store,
 * @author Lars Simon Zehnder
 * ---------------------------------------------------------
 **/
template <typename Super, typename PostOutType>
void POST <Super, PostOutType>::Output::store(const unsigned int& m,
                                              Node& node)
{
   Super::Output::store(m, node);
   if (m >= node.BURNIN)
   {
      const unsigned int index = m - node.BURNIN;
      post.store(index, node.hyperPar);
   }
}

// ========================================================
// POST mixin layer
// --------------------------------------------------------

/**
 * --------------------------------------------------------
 * POST<Super>::POST
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
 *          mixins. Calls the constructor of its Super
 *          layer in initializing list.
 * @see Super, FinmixData, FinmixModel, FinmixPrior,
 *      FinmixMCMC, Rcpp::S4
 * @author  Lars Simon Zehnder
 * -------------------------------------------------------
 **/
template <typename Super, typename PostOutType>
POST <Super, PostOutType>::POST (const FinmixData& data,
                                 const FinmixModel& model, const FinmixPrior& prior,
                                 const FinmixMCMC& mcmc, Rcpp::S4& classS4) :
   Super(data, model, prior, mcmc, classS4),
   node(data, model, prior, mcmc),
   output(classS4)
{
}

/**
 * -------------------------------------------------------
 * POST<Super>::update
 * @brief   Triggers the update process for each step
 *          the sampler. Passes responsibility to 'Node's
 *          'update()' method.
 * @see Node::update
 * @author  Lars Simon Zehnder
 * -------------------------------------------------------
 **/
template <typename Super, typename PostOutType>
void POST <Super, PostOutType>::update()
{
   node.update();
}

/**
 * -------------------------------------------------------
 * POST<Super>::store
 * @brief   Triggers the store process for each step of
 *          the sampler. Passes responsibility to 'Output's
 *          'store()' method.
 * @see Output::store
 * @author  Lars Simon Zehnder
 * -------------------------------------------------------
 **/
template <typename Super, typename PostOutType>
void POST <Super, PostOutType>::store(const unsigned int& m)
{
   output.store(m, node);
}
#endif
