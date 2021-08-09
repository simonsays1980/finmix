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
#ifndef IND_H
#define IND_H

#include <RcppArmadillo.h>

// ===================================================================
// IND mixin layer
// -------------------------------------------------------------------
/*
 * @brief   Mixin layer to implement the collaboration between 'Node'
 *          and 'Output' objects in case of Gibbs sampling with 
 *          indicators.
 * @par Super   next mixin layer in the application
 * @detail  Any implemented mixin layer describes the whole collabo-
 *          ration between 'Node' and 'Output' object to perform a
 *          Gibbs sampling of posterior parameters. The IND mixin
 *          layer refines thereby its 'Super' class by defining 
 *          new inner mixins 'Node' and 'Output' with additional
 *          variables needed to perform all actions for Gibbs 
 *          sampling with indicators. These are e.g. variables 
 *          neded to perform the permutations for random permutation
 *          Gibbs sampling. 
 * @see FIX, HIER, POST, ADAPTER, BASE
 * @author Lars SImon Zehnder
 * ------------------------------------------------------------------
 */
template <typename Super> 
class IND : public Super {
	public:
        /**
         * ---------------------------------------------------------
         * Node mixin 
         * ---------------------------------------------------------
         * 
         * @brief   Holds all variables and method to perform the 
         *          steps of a Gibbs sampler. 
         * @detail  This class inherits directly from the 'Super'
         *          class' 'Node' mixin. It defines the new var-
         *          iable 'swapIndex' needed for random permutation
         *          Gibbs sampling including classification sampling
         *          and permutating. The workhorse of this mixin is 
         *          the inherited virtual method 'update()' that 
         *          performs the update step and calls any 'update()'
         *          function of related classes.
         * @see FIX, HIER, POST, ADAPTER, BASE
         * --------------------------------------------------------
         */
		class Node : public Super::Node {
			public: 
				arma::urowvec swapIndex;

				Node (const FinmixData&,
					const FinmixModel&,
					const FinmixPrior&,
					const FinmixMCMC&);
				virtual ~Node () {}
				virtual void update (); 
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
         *          tional information from sampling indicators.
         *          Reusable functionality is inherited from 
         *          'Super's 'Output' class. The workhorse of this
         *          inner mixin is the 'store()' method that per-
         *          forms the storing process thereby calling all
         *          'store()' methods of related classes.
         * @see FIX, BASE, HIER, POST, ADAPTER
         * ------------------------------------------------------
         */
		class Output : public Super::Output {
			public:
				arma::mat* weight;
				arma::vec* cdpost;
				arma::vec* entropy;
				arma::ivec* ST;
				arma::imat* S;
				arma::imat* NK;
				arma::ivec* clust;
				Output (Rcpp::S4&);
				virtual ~Output () {}
				virtual void store (const unsigned int&,
					Node&);
		};
		Node node;
		Output output;

		IND (const FinmixData&, const FinmixModel&,
			const FinmixPrior&, const FinmixMCMC&,
			Rcpp::S4&);
		virtual ~IND () {}
		virtual void update ();
		virtual void store (const unsigned int&);
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
template <typename Super>
IND <Super>::Node::Node (const FinmixData& data, 
	const FinmixModel& model, const FinmixPrior& prior,
	const FinmixMCMC& mcmc) :
		Super::Node(data, model, prior, mcmc),
		swapIndex(model.K) {}

/**
 * -----------------------------------------------------------
 * Node::update
 * @brief   Updates the 'node' object.
 * @detail  Virtual. Calls 'Super::Node::update()', to perform
 *          any updates on parameters and then starts random 
 *          permutation of sampled parameters and indicators.
 * @see Super::Node::update()
 * @author  Lars Simon Zehnder
 * -----------------------------------------------------------
 **/
template <typename Super>
void IND <Super>::Node::update () 
{
	Super::Node::update();
	if (Super::Node::RANPERM && arma::sum(Super::Node::compIndex2) !=
		Super::Node::K) {
		Super::Node::par.weight(Super::Node::compIndex) = 
			Super::Node::par.weight(Super::Node::permIndex);
   		swapIndex = arma::sort_index(Super::Node::permIndex).t();
        for(unsigned int i = 0; i < Super::Node::N; ++i) {
			Super::Node::S(i) = (int) swapIndex((unsigned int) 
				(Super::Node::S(i) - 1)) + 1;
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
 *          memory and fix the size.
 * @see FIX::Output::Output, HIER::Output::Output,
 *      Rcpp::S4, ?S4 (in R), arma::mat
 * @author  Lars Simon Zehnder
 * ----------------------------------------------------------
 **/
template <typename Super>
IND <Super>::Output::Output (Rcpp::S4& classS4) : 
	Super::Output(classS4) 
{
	Rcpp::NumericMatrix tmpWeight((SEXP) classS4.slot("weight"));
	Rcpp::List tmpList((SEXP) classS4.slot("log"));
	Rcpp::NumericVector tmpCDPost((SEXP) tmpList["cdpost"]);
	Rcpp::NumericVector tmpEntropy((SEXP) classS4.slot("entropy"));
	Rcpp::IntegerVector tmpST((SEXP) classS4.slot("ST"));
	Rcpp::IntegerMatrix tmpS((SEXP) classS4.slot("S"));
	Rcpp::IntegerMatrix tmpNK((SEXP) classS4.slot("NK"));
	Rcpp::IntegerVector tmpClust((SEXP) classS4.slot("clust"));
	const unsigned int tmpM = tmpWeight.nrow();
	const unsigned int K = tmpWeight.ncol();
	const unsigned int N = tmpS.nrow();
	const unsigned int STORES = tmpS.ncol();
	weight = new arma::mat(tmpWeight.begin(), tmpM, K, false, true);
	cdpost = new arma::vec(tmpCDPost.begin(), tmpM, false, true);
	entropy = new arma::vec(tmpEntropy.begin(), tmpM, false, true);
	ST = new arma::ivec(tmpST.begin(), tmpM, false, true);
	S = new arma::imat(tmpS.begin(), N, STORES, false, true);
	NK = new arma::imat(tmpNK.begin(), tmpM, K, false, true);
	clust = new arma::ivec(tmpClust.begin(), N, false, true);
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
 *          phase or the sampling phase, if indicators
 *          should be stored at all, etc.
 * @see FIX::Output::store, HIER::Output::store
 * @author Lars Simon Zehnder
 * ---------------------------------------------------------
 **/
template <typename Super> 
void IND <Super>::Output::store (const unsigned int& m,
	Node& node)
{
	Super::Output::store(m,node);
	if (m >= node.BURNIN) {
		const unsigned int index = m - node.BURNIN;
		(*weight).row(index) = node.par.weight;
		(*cdpost)(index) = node.log.cdpost;
		(*entropy)(index) = node.log.entropy;
		(*ST)(index) = node.S(node.N - 1);
		if(index >= node.M - node.STORES) {
			if (!node.STARTPAR && index != node.M - 1) {
				(*S).col(index - (node.M - node.STORES) + 1) = node.S;
			}
			if (node.STARTPAR){
				(*S).col(index - (node.M - node.STORES)) = node.S;
			}
		}
		(*NK).row(index) = arma::conv_to<arma::irowvec>::from
			(node.hyperPar.weightPost - node.hyperPar.weightStart);
		if (m == node.BURNIN) {
			node.log.maxcdpost = node.log.cdpost - 1;		
		}
		if (node.log.cdpost > node.log.maxcdpost) {
			(*clust) = node.S;
		}
	}
}

// ========================================================
// IND mixin layer
// --------------------------------------------------------

/**
 * --------------------------------------------------------
 * IND<Super>::IND
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
template <typename Super>
IND <Super>::IND (const FinmixData& data, const FinmixModel& model,
	const FinmixPrior& prior, const FinmixMCMC& mcmc,
	Rcpp::S4& classS4) :
		Super(data, model, prior, mcmc, classS4),
		node(data, model, prior, mcmc),
		output(classS4) {}

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
template <typename Super>
void IND <Super>::update ()
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
template <typename Super>
void IND <Super>::store (const unsigned int& m)
{
	output.store(m, node);
}  
#endif
