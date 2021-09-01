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

// [[Rcpp::export]]

Rcpp::NumericMatrix swap_cc(Rcpp::NumericMatrix values, Rcpp::IntegerMatrix index) {
    /* If dimensions of both arguments do not agree throw an exception */
    if(values.nrow() != index.nrow() || values.ncol() != index.ncol()) {
        throw Rcpp::exception("Matrix dimensions disagree.");
    }
    /* Do not reuse memory from R as otherwise existing objects 
     * get manipulated */
    const unsigned int K = values.ncol();
    const unsigned int M = values.nrow();
    arma::mat values_arma(values.begin(), M, K, true, true);
    arma::imat index_arma(index.begin(), M, K, true, true);
    arma::mat values_copy(M, K);
    arma::umat index_umat = arma::conv_to<arma::umat>::from(index_arma) - 1;
    arma::uvec row_index(1);
    arma::urowvec swap_index(K);
    for(unsigned int i = 0; i < M; ++i) {
        row_index.at(0) = i;
        swap_index = index_umat.row(i);
        values_copy.row(i) = 
            values_arma.submat(row_index, swap_index);
    }
    return Rcpp::wrap(values_copy);
}

// [[Rcpp::export]]
Rcpp::NumericVector swap_3d_cc(Rcpp::NumericVector values, Rcpp::IntegerMatrix index)
{
    Rcpp::IntegerVector valDim = values.attr("dim");
    const unsigned int M = valDim[0];
    const unsigned int r = valDim[1];
    const unsigned int K = valDim[2];
    /* If dimensions of both arguments do not agree thrw an exception */
    if ( M != (unsigned)index.nrow() || K != (unsigned)index.ncol()) {
        throw Rcpp::exception("Matrix dimensions disagree.");
    }
    arma::cube values_arma(values.begin(), M, r, K, false, true);
    arma::imat index_arma(index.begin(), M, K, false, true);
    arma::cube output(M, r, K);
    output.fill(0.0);
    arma::umat index_umat = arma::conv_to<arma::umat>::from(index_arma) - 1;
    arma::umat ik(M, K);
    arma::ucube ikr(M, 1, K);
    arma::ucube ikr2(M, r, K);
    arma::cube ikr3(M, r, K);
    for (unsigned int k = 0; k < K; ++k) {
        ik = (index_arma - 1) == k;
        ikr.slices(0, K - 1) = ik;
        ikr2 = arma::resize(ikr, M, r, K);
        for (unsigned int rr = 1; rr < r; ++rr) {
            ikr2.tube(0, rr, M - 1, rr) = ikr2.tube(0, 0, M - 1, 0);
        }
        ikr3 = arma::conv_to<arma::cube>::from(ikr2);
        ikr3 %= values_arma;
        for (unsigned int l = 0; l < K; ++l) {
            output.slice(k) += ikr3.slice(l);    
        }
        ik.fill(0);
        ikr.fill(0);
        ikr2.fill(0);
        ikr3.fill(0);
    }
    return Rcpp::wrap(output);
}

// [[Rcpp::export]]

Rcpp::IntegerMatrix swapInteger_cc(Rcpp::IntegerMatrix values, Rcpp::IntegerMatrix index) {
    /* If dimensions of both arguments do not agree throw an exception */
    if(values.nrow() != index.nrow() || values.ncol() != index.ncol()) {
        throw Rcpp::exception("Matrix dimensions disagree.");
    }
    /* Do not reuse memory from R as otherwise existing objects 
     * get manipulated */
    const unsigned int K = values.ncol();
    const unsigned int M = values.nrow();
    arma::imat values_arma(values.begin(), M, K, true, true);
    arma::imat index_arma(index.begin(), M, K, true, true);
    arma::imat values_copy(M, K);
    arma::umat index_umat = arma::conv_to<arma::umat>::from(index_arma) - 1;
    arma::uvec row_index(1);
    arma::urowvec swap_index(K);
    for(unsigned int i = 0; i < M; ++i) {
        row_index.at(0) = i;
        swap_index = index_umat.row(i);
        values_copy.row(i) = 
            values_arma.submat(row_index, swap_index);
    }
    return Rcpp::wrap(values_copy);
}

// [[Rcpp::export]]

Rcpp::IntegerMatrix swapInd_cc(Rcpp::IntegerMatrix values, Rcpp::IntegerMatrix index) {
    /* If dimensions of both arguments do not agree throw an exception */
    if (values.ncol() != index.nrow()) {
        throw Rcpp::exception("Matrix dimensions disagree.");
    }
    /* Reuse memory from R */
    const unsigned int N = values.nrow();
    const unsigned int STORES = values.ncol();
    const unsigned int M = index.nrow();
    const unsigned int K = index.ncol();
    arma::imat values_arma(values.begin(), N, STORES, true, true);
    arma::imat index_arma(index.begin(), M, K, true, true);
    arma::imat values_copy(N, STORES);
    for (unsigned int s = 0; s < STORES; ++s) {
        for (unsigned int i = 0; i < N; ++i) {
            values_copy(i, s) = (int) index_arma(s, (unsigned int)
                    values_arma(i, s) - 1);
        }
    }
    return Rcpp::wrap(values_copy);
}

// [[Rcpp::export]]

Rcpp::IntegerVector swapST_cc(Rcpp::IntegerVector values, Rcpp::IntegerMatrix index) {
    /* If dimensions of both arguments do not agree throw an exception */
    if (values.size() != index.nrow()) {
        throw Rcpp::exception("Matrix dimensions disagree.");
    }
    /* Reuse memory from R */
    const unsigned int M = values.size();
    const unsigned int K = index.ncol();
    arma::ivec values_arma(values.begin(), M, false, true);
    arma::imat index_arma(index.begin(), M, K, false, true);
    arma::ivec values_copy(M);
    for(unsigned int i = 0; i < M; ++i) {
        values_copy(i) = index_arma(i, (unsigned int)
                values_arma(i) - 1);
    }
    return Rcpp::wrap(values_copy);
}

// [[Rcpp::export]]

Rcpp::NumericMatrix ldgamma_cc(Rcpp::NumericMatrix values, Rcpp::NumericVector shape,
        Rcpp::NumericVector rate)
{
    /* Reuse memory from R */
    const unsigned int M    = values.nrow();
    const unsigned int K    = values.ncol();
    arma::mat arma_values(values.begin(), M, K, false, true);
    arma::vec arma_shape(shape.begin(), K, false, true);
    arma::vec arma_rate(rate.begin(), K, false, true);
    arma::mat arma_return(M, K);
    arma_return = ldgamma(arma_values, arma_shape, arma_rate);
    return Rcpp::wrap(arma_return);
}

// [[Rcpp::export]]

arma::mat dgamma_cc(Rcpp::NumericMatrix values, Rcpp::NumericVector shape,
        Rcpp::NumericVector rate)
{
    /* Reuse memory from R */
    const unsigned int M    = values.nrow();
    const unsigned int K    = values.ncol();
    arma::mat arma_values(values.begin(), M, K, false, true);
    arma::vec arma_shape(shape.begin(), K, false, true);
    arma::vec arma_rate(rate.begin(), K, false, true);
    arma::mat arma_return(M, K);
    arma_return = exp(ldgamma(arma_values, arma_shape, arma_rate));
    return arma_return;
}

// [[Rcpp::export]]

Rcpp::NumericVector lddirichlet_cc(Rcpp::NumericMatrix values, Rcpp::NumericVector par)
{
    /* Reuse memory from R */
    const unsigned int M    = values.nrow();
    const unsigned int K    = values.ncol();
    arma::mat arma_values(values.begin(), M, K, false, true);
    arma::vec arma_par(par.begin(), K, false, true);
    arma::vec arma_return(M);
    arma_return = lddirichlet(arma_values, arma_par);
    return Rcpp::wrap(arma_return);
}

// [[Rcpp::export]]

arma::vec ddirichlet_cc(Rcpp::NumericMatrix values, Rcpp::NumericVector par)
{
    /* Reuse memory from R */
    const unsigned int M    = values.nrow();
    const unsigned int K    = values.ncol();
    arma::mat arma_values(values.begin(), M, K, false, true);
    arma::vec arma_par(par.begin(), K, false, true);
    arma::vec arma_return(M);
    arma_return = arma::exp(lddirichlet(arma_values, arma_par));
    return arma_return;
}

// [[Rcpp::export]] 

arma::imat hungarian_cc(const arma::mat cost) 
{
    arma::umat indM = hungarian(cost);
    return arma::conv_to<arma::imat>::from(indM);
}

// [[Rcpp::export]]

Rcpp::List moments_cc(Rcpp::S4 classS4) 
{
    Rcpp::S4 model      = Rcpp::as<Rcpp::S4>((SEXP) classS4.slot("model"));
    const bool indicfix = Rcpp::as<bool>((SEXP) model.slot("indicfix"));
    if (indicfix) {
        return Rcpp::wrap(moments_fix_cc(classS4));
    } else {
        return Rcpp::wrap(moments_ind_cc(classS4));
    }
}

// [[Rcpp::export]]

Rcpp::List permmoments_cc(Rcpp::S4 classS4) 
{
    Rcpp::S4 model      = Rcpp::as<Rcpp::S4>((SEXP) classS4.slot("model"));
    const bool indicfix = Rcpp::as<bool>((SEXP) model.slot("indicfix"));
    if (indicfix) {
        return Rcpp::wrap(permmoments_fix_cc(classS4));
    } else {
        return Rcpp::wrap(permmoments_ind_cc(classS4));
    }
}

