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

#ifndef __FINMIX_MOMENTS_H__
#define __FINMIX_MOMENTS_H__

#include <RcppArmadillo.h> 

Rcpp::List moments_fix_cc (Rcpp::S4 classS4) 
{
    Rcpp::List parList              = Rcpp::as<Rcpp::List>((SEXP) classS4.slot("par"));
    Rcpp::NumericVector tmpMu       = Rcpp::as<Rcpp::NumericVector>((SEXP) parList["mu"]);
    Rcpp::NumericVector tmpSigma    = Rcpp::as<Rcpp::NumericVector>((SEXP) parList["sigma"]);
    Rcpp::IntegerVector tmpMuDim    = tmpMu.attr("dim");
    Rcpp::IntegerVector tmpSigmaDim = tmpSigma.attr("dim");
    const unsigned int M            = tmpMuDim[0];
    const unsigned int r            = tmpMuDim[1];
    const unsigned int K            = tmpMuDim[2];
    const unsigned int s            = tmpSigmaDim[1];
    const unsigned int ij           = R::choose(r, 2);
    arma::cube mu                   = arma::cube(tmpMu.begin(), M, r, K, false, true);
    arma::cube sigma                = arma::cube(tmpSigma.begin(), M, s, K, false, true);
    Rcpp::S4 model                  = Rcpp::as<Rcpp::S4>((SEXP) classS4.slot("model"));
    arma::vec weight                = Rcpp::as<arma::vec>((SEXP) model.slot("weight"));   
    arma::vec means(r);
    arma::vec tmp(K);
    arma::mat var(r, r);
    arma::mat W(r, r);
    arma::mat B(r, r);
    arma::mat cd(r, r);
    arma::mat corr(r, r);
    arma::rowvec tmp2;
    arma::vec d;
    double Rtr                      = 0.0;
    double Rdet                     = 0.0;
    arma::vec zm(4);
    arma::mat higher(r, 4);
    arma::vec sigmavec(K);
    arma::vec cm(K);
    arma::vec skewness(r);
    arma::vec kurtosis(r);

    // Output containers 
    arma::mat meanOut(M, r);
    arma::vec RtrOut(M);
    arma::vec RdetOut(M);
    arma::mat corrOut(M, ij);
    arma::mat varOut(M, r);
    arma::mat skewnessOut(M, r);
    arma::mat kurtosisOut(M, r);

    // Permutation matrix
    arma::umat perm(ij, 2);
    unsigned int index = 0;
    do {
        for (unsigned int i = 0; i < r - 1; ++i) {
            for (unsigned int j = i + 1; j < r; ++j) {
                perm(index, 0) = i;
                perm(index, 1) = j;
                ++index;
            }
        }
    } while (index < ij - 1);

    for (unsigned int i = 0; i < M; ++i) {       
        higher.fill(0.0);
        for (unsigned int rr = 0; rr < r; ++rr) {
            tmp = mu.tube(i, rr, i, rr);
            means(rr) = arma::as_scalar(weight.t() * tmp);            
        }
        var.fill(0.0);
        W.fill(0.0);
        B.fill(0.0);
        cd.fill(0.0);
        corr.fill(0.0);
        
        for (unsigned int k = 0; k < K; ++k) {
            tmp2 = mu.slice(k)(arma::span(i), arma::span());
            var = var + tmp2.t() * tmp2 
                + qinmatr(sigma.slice(k)(arma::span(i),arma::span())) 
                * weight(k);
            W   = W + qinmatr(sigma.slice(k)(arma::span(i),arma::span()))
                * weight(k);
            d   = arma::trans(mu.slice(k)(arma::span(i), arma::span()))
                - means;
            B   = B + d * d.t() * weight(k);
        }
        var     = var - means * means.t();
        cd      = arma::diagmat(1.0 / arma::sqrt(arma::diagvec(var)));
        corr    = cd * var * cd;
        Rtr     = 1 - arma::trace(W) / arma::trace(var);
        Rdet    = 1 - std::log(arma::det(W)) / std::log(arma::det(var));
        zm.fill(0.0);        
        zm(1)   = 1.0;
        zm(3)   = std::exp(std::log(1.0) + std::log(3.0));        
        for (unsigned int m = 0; m < 4; ++m) {
            for (unsigned int rr = 0; rr < r; ++rr) {
                for (unsigned int k = 0; k < K; ++k) {
                    sigmavec(k) = qinmatr(sigma.slice(k)(arma::span(i), arma::span()))(rr,rr);                    
                }
                tmp = mu.tube(i, rr, i, rr) - means(rr); 
                higher(rr, m)   = arma::as_scalar(weight.t() * arma::pow(tmp, m + 1));
                for (unsigned int n = 0; n < (m + 1); ++n) {
                    arma::vec ss = arma::pow(tmp, m - n);
                    cm = arma::pow(tmp, m - n) % arma::pow(sigmavec, (n + 1) / 2)
                        * zm[n];
                    higher(rr, m)   = higher(rr, m) + R::choose(m + 1, n + 1) 
                        * arma::as_scalar(weight.t() * cm);                    
                }
            }
        }
        skewness    = higher.col(2) / arma::pow(higher.col(1), 1.5);
        kurtosis    = higher.col(3) / arma::pow(higher.col(1), 2);
        meanOut.row(i) = arma::trans(means);
        RtrOut(i)   = Rtr;
        RdetOut(i)  = Rdet;
        for (unsigned int j = 0; j < ij; ++j) {
            corrOut(i, j) = corr(perm(j, 0), perm(j, 1));
        }
        for (unsigned int rr = 0; rr < r; ++rr) {
            varOut(i, rr) = var(rr, rr);
        }
        skewnessOut.row(i) = arma::trans(skewness);
        kurtosisOut.row(i) = arma::trans(kurtosis);
    }

    return Rcpp::List::create(Rcpp::Named("Rtr", RtrOut),
            Rcpp::Named("mean", meanOut),
            Rcpp::Named("Rdet", RdetOut),
            Rcpp::Named("corr", corrOut),
            Rcpp::Named("var", varOut),
            Rcpp::Named("skewness", skewnessOut),
            Rcpp::Named("kurtosis", kurtosisOut));
}

Rcpp::List moments_ind_cc (Rcpp::S4 classS4) 
{
    Rcpp::List parList              = Rcpp::as<Rcpp::List>((SEXP) classS4.slot("par"));
    Rcpp::NumericVector tmpMu       = Rcpp::as<Rcpp::NumericVector>((SEXP) parList["mu"]);
    Rcpp::NumericVector tmpSigma    = Rcpp::as<Rcpp::NumericVector>((SEXP) parList["sigma"]);
    Rcpp::IntegerVector tmpMuDim    = tmpMu.attr("dim");
    Rcpp::IntegerVector tmpSigmaDim = tmpSigma.attr("dim");
    const unsigned int M            = tmpMuDim[0];
    const unsigned int r            = tmpMuDim[1];
    const unsigned int K            = tmpMuDim[2];
    const unsigned int s            = tmpSigmaDim[1];
    const unsigned int ij           = R::choose(r, 2);
    arma::cube mu                   = arma::cube(tmpMu.begin(), M, r, K, false, true);
    arma::cube sigma                = arma::cube(tmpSigma.begin(), M, s, K, false, true);
    arma::mat weights               = Rcpp::as<arma::mat>((SEXP) classS4.slot("weight"));  
    arma::vec weight(K);
    arma::vec means(r);
    arma::vec tmp(K);
    arma::mat var(r, r);
    arma::mat W(r, r);
    arma::mat B(r, r);
    arma::mat cd(r, r);
    arma::mat corr(r, r);
    arma::rowvec tmp2;
    arma::vec d;
    double Rtr                      = 0.0;
    double Rdet                     = 0.0;
    arma::vec zm(4);
    arma::mat higher(r, 4);
    arma::vec sigmavec(K);
    arma::vec cm(K);
    arma::vec skewness(r);
    arma::vec kurtosis(r);

    // Output containers 
    arma::mat meanOut(M, r);
    arma::vec RtrOut(M);
    arma::vec RdetOut(M);
    arma::mat corrOut(M, ij);
    arma::mat varOut(M, r);
    arma::mat skewnessOut(M, r);
    arma::mat kurtosisOut(M, r);

    // Permutation matrix
    arma::umat perm(ij, 2);
    unsigned int index = 0;
    do {
        for (unsigned int i = 0; i < r - 1; ++i) {
            for (unsigned int j = i + 1; j < r; ++j) {
                perm(index, 0) = i;
                perm(index, 1) = j;
                ++index;
            }
        }
    } while (index < ij - 1);

    for (unsigned int i = 0; i < M; ++i) {       
        higher.fill(0.0);
        weight = arma::trans(weights.row(i));
        for (unsigned int rr = 0; rr < r; ++rr) {
            tmp = mu.tube(i, rr, i, rr);
            means(rr) = arma::as_scalar(weight.t() * tmp);            
        }
        var.fill(0.0);
        W.fill(0.0);
        B.fill(0.0);
        cd.fill(0.0);
        corr.fill(0.0);
        
        for (unsigned int k = 0; k < K; ++k) {
            tmp2 = mu.slice(k)(arma::span(i), arma::span());
            var = var + tmp2.t() * tmp2 
                + qinmatr(sigma.slice(k)(arma::span(i),arma::span())) 
                * weight(k);
            W   = W + qinmatr(sigma.slice(k)(arma::span(i),arma::span()))
                * weight(k);
            d   = arma::trans(mu.slice(k)(arma::span(i), arma::span()))
                - means;
            B   = B + d * d.t() * weight(k);
        }
        var     = var - means * means.t();
        cd      = arma::diagmat(1.0 / arma::sqrt(arma::diagvec(var)));
        corr    = cd * var * cd;
        Rtr     = 1 - arma::trace(W) / arma::trace(var);
        Rdet    = 1 - std::log(arma::det(W)) / std::log(arma::det(var));
        zm.fill(0.0);        
        zm(1)   = 1.0;
        zm(3)   = std::exp(std::log(1.0) + std::log(3.0));        
        for (unsigned int m = 0; m < 4; ++m) {
            for (unsigned int rr = 0; rr < r; ++rr) {
                for (unsigned int k = 0; k < K; ++k) {
                    sigmavec(k) = qinmatr(sigma.slice(k)(arma::span(i), arma::span()))(rr,rr);                    
                }
                tmp = mu.tube(i, rr, i, rr) - means(rr); 
                higher(rr, m)   = arma::as_scalar(weight.t() * arma::pow(tmp, m + 1));
                for (unsigned int n = 0; n < (m + 1); ++n) {
                    arma::vec ss = arma::pow(tmp, m - n);
                    cm = arma::pow(tmp, m - n) % arma::pow(sigmavec, (n + 1) / 2)
                        * zm[n];
                    higher(rr, m)   = higher(rr, m) + R::choose(m + 1, n + 1) 
                        * arma::as_scalar(weight.t() * cm);                    
                }
            }
        }
        skewness    = higher.col(2) / arma::pow(higher.col(1), 1.5);
        kurtosis    = higher.col(3) / arma::pow(higher.col(1), 2);
        meanOut.row(i) = arma::trans(means);
        RtrOut(i)   = Rtr;
        RdetOut(i)  = Rdet;
        for (unsigned int j = 0; j < ij; ++j) {
            corrOut(i, j) = corr(perm(j, 0), perm(j, 1));
        }
        for (unsigned int rr = 0; rr < r; ++rr) {
            varOut(i, rr) = var(rr, rr);
        }
        skewnessOut.row(i) = arma::trans(skewness);
        kurtosisOut.row(i) = arma::trans(kurtosis);
    }

    return Rcpp::List::create(Rcpp::Named("Rtr", RtrOut),
            Rcpp::Named("mean", meanOut),
            Rcpp::Named("Rdet", RdetOut),
            Rcpp::Named("corr", corrOut),
            Rcpp::Named("var", varOut),
            Rcpp::Named("skewness", skewnessOut),
            Rcpp::Named("kurtosis", kurtosisOut));
}

Rcpp::List permmoments_fix_cc (Rcpp::S4 classS4) 
{
    Rcpp::List parList              = Rcpp::as<Rcpp::List>((SEXP) classS4.slot("parperm"));
    Rcpp::NumericVector tmpMu       = Rcpp::as<Rcpp::NumericVector>((SEXP) parList["mu"]);
    Rcpp::NumericVector tmpSigma    = Rcpp::as<Rcpp::NumericVector>((SEXP) parList["sigma"]);
    Rcpp::IntegerVector tmpMuDim    = tmpMu.attr("dim");
    Rcpp::IntegerVector tmpSigmaDim = tmpSigma.attr("dim");
    const unsigned int M            = tmpMuDim[0];
    const unsigned int r            = tmpMuDim[1];
    const unsigned int K            = tmpMuDim[2];
    const unsigned int s            = tmpSigmaDim[1];
    const unsigned int ij           = R::choose(r, 2);
    arma::cube mu                   = arma::cube(tmpMu.begin(), M, r, K, false, true);
    arma::cube sigma                = arma::cube(tmpSigma.begin(), M, s, K, false, true);
    Rcpp::S4 model                  = Rcpp::as<Rcpp::S4>((SEXP) classS4.slot("model"));
    arma::vec weight                = Rcpp::as<arma::vec>((SEXP) model.slot("weight"));   
    arma::vec means(r);
    arma::vec tmp(K);
    arma::mat var(r, r);
    arma::mat W(r, r);
    arma::mat B(r, r);
    arma::mat cd(r, r);
    arma::mat corr(r, r);
    arma::rowvec tmp2;
    arma::vec d;
    double Rtr                      = 0.0;
    double Rdet                     = 0.0;
    arma::vec zm(4);
    arma::mat higher(r, 4);
    arma::vec sigmavec(K);
    arma::vec cm(K);
    arma::vec skewness(r);
    arma::vec kurtosis(r);

    // Output containers 
    arma::mat meanOut(M, r);
    arma::vec RtrOut(M);
    arma::vec RdetOut(M);
    arma::mat corrOut(M, ij);
    arma::mat varOut(M, r);
    arma::mat skewnessOut(M, r);
    arma::mat kurtosisOut(M, r);

    // Permutation matrix
    arma::umat perm(ij, 2);
    unsigned int index = 0;
    do {
        for (unsigned int i = 0; i < r - 1; ++i) {
            for (unsigned int j = i + 1; j < r; ++j) {
                perm(index, 0) = i;
                perm(index, 1) = j;
                ++index;
            }
        }
    } while (index < ij - 1);

    for (unsigned int i = 0; i < M; ++i) {       
        higher.fill(0.0);
        for (unsigned int rr = 0; rr < r; ++rr) {
            tmp = mu.tube(i, rr, i, rr);
            means(rr) = arma::as_scalar(weight.t() * tmp);            
        }
        var.fill(0.0);
        W.fill(0.0);
        B.fill(0.0);
        cd.fill(0.0);
        corr.fill(0.0);
        
        for (unsigned int k = 0; k < K; ++k) {
            tmp2 = mu.slice(k)(arma::span(i), arma::span());
            var = var + tmp2.t() * tmp2 
                + qinmatr(sigma.slice(k)(arma::span(i),arma::span())) 
                * weight(k);
            W   = W + qinmatr(sigma.slice(k)(arma::span(i),arma::span()))
                * weight(k);
            d   = arma::trans(mu.slice(k)(arma::span(i), arma::span()))
                - means;
            B   = B + d * d.t() * weight(k);
        }
        var     = var - means * means.t();
        cd      = arma::diagmat(1.0 / arma::sqrt(arma::diagvec(var)));
        corr    = cd * var * cd;
        Rtr     = 1 - arma::trace(W) / arma::trace(var);
        Rdet    = 1 - std::log(arma::det(W)) / std::log(arma::det(var));
        zm.fill(0.0);        
        zm(1)   = 1.0;
        zm(3)   = std::exp(std::log(1.0) + std::log(3.0));        
        for (unsigned int m = 0; m < 4; ++m) {
            for (unsigned int rr = 0; rr < r; ++rr) {
                for (unsigned int k = 0; k < K; ++k) {
                    sigmavec(k) = qinmatr(sigma.slice(k)(arma::span(i), arma::span()))(rr,rr);                    
                }
                tmp = mu.tube(i, rr, i, rr) - means(rr); 
                higher(rr, m)   = arma::as_scalar(weight.t() * arma::pow(tmp, m + 1));
                for (unsigned int n = 0; n < (m + 1); ++n) {
                    arma::vec ss = arma::pow(tmp, m - n);
                    cm = arma::pow(tmp, m - n) % arma::pow(sigmavec, (n + 1) / 2)
                        * zm[n];
                    higher(rr, m)   = higher(rr, m) + R::choose(m + 1, n + 1) 
                        * arma::as_scalar(weight.t() * cm);                    
                }
            }
        }
        skewness    = higher.col(2) / arma::pow(higher.col(1), 1.5);
        kurtosis    = higher.col(3) / arma::pow(higher.col(1), 2);
        meanOut.row(i) = arma::trans(means);
        RtrOut(i)   = Rtr;
        RdetOut(i)  = Rdet;
        for (unsigned int j = 0; j < ij; ++j) {
            corrOut(i, j) = corr(perm(j, 0), perm(j, 1));
        }
        for (unsigned int rr = 0; rr < r; ++rr) {
            varOut(i, rr) = var(rr, rr);
        }
        skewnessOut.row(i) = arma::trans(skewness);
        kurtosisOut.row(i) = arma::trans(kurtosis);
    }

    return Rcpp::List::create(Rcpp::Named("Rtr", RtrOut),
            Rcpp::Named("mean", meanOut),
            Rcpp::Named("Rdet", RdetOut),
            Rcpp::Named("corr", corrOut),
            Rcpp::Named("var", varOut),
            Rcpp::Named("skewness", skewnessOut),
            Rcpp::Named("kurtosis", kurtosisOut));
}

Rcpp::List permmoments_ind_cc (Rcpp::S4 classS4) 
{
    Rcpp::List parList              = Rcpp::as<Rcpp::List>((SEXP) classS4.slot("parperm"));
    Rcpp::NumericVector tmpMu       = Rcpp::as<Rcpp::NumericVector>((SEXP) parList["mu"]);
    Rcpp::NumericVector tmpSigma    = Rcpp::as<Rcpp::NumericVector>((SEXP) parList["sigma"]);
    Rcpp::IntegerVector tmpMuDim    = tmpMu.attr("dim");
    Rcpp::IntegerVector tmpSigmaDim = tmpSigma.attr("dim");
    const unsigned int M            = tmpMuDim[0];
    const unsigned int r            = tmpMuDim[1];
    const unsigned int K            = tmpMuDim[2];
    const unsigned int s            = tmpSigmaDim[1];
    const unsigned int ij           = R::choose(r, 2);
    arma::cube mu                   = arma::cube(tmpMu.begin(), M, r, K, false, true);
    arma::cube sigma                = arma::cube(tmpSigma.begin(), M, s, K, false, true);
    arma::mat weights               = Rcpp::as<arma::mat>((SEXP) classS4.slot("weightperm"));  
    arma::vec weight(K);
    arma::vec means(r);
    arma::vec tmp(K);
    arma::mat var(r, r);
    arma::mat W(r, r);
    arma::mat B(r, r);
    arma::mat cd(r, r);
    arma::mat corr(r, r);
    arma::rowvec tmp2;
    arma::vec d;
    double Rtr                      = 0.0;
    double Rdet                     = 0.0;
    arma::vec zm(4);
    arma::mat higher(r, 4);
    arma::vec sigmavec(K);
    arma::vec cm(K);
    arma::vec skewness(r);
    arma::vec kurtosis(r);

    // Output containers 
    arma::mat meanOut(M, r);
    arma::vec RtrOut(M);
    arma::vec RdetOut(M);
    arma::mat corrOut(M, ij);
    arma::mat varOut(M, r);
    arma::mat skewnessOut(M, r);
    arma::mat kurtosisOut(M, r);

    // Permutation matrix
    arma::umat perm(ij, 2);
    unsigned int index = 0;
    do {
        for (unsigned int i = 0; i < r - 1; ++i) {
            for (unsigned int j = i + 1; j < r; ++j) {
                perm(index, 0) = i;
                perm(index, 1) = j;
                ++index;
            }
        }
    } while (index < ij - 1);

    for (unsigned int i = 0; i < M; ++i) {       
        higher.fill(0.0);
        weight = arma::trans(weights.row(i));
        for (unsigned int rr = 0; rr < r; ++rr) {
            tmp = mu.tube(i, rr, i, rr);
            means(rr) = arma::as_scalar(weight.t() * tmp);            
        }
        var.fill(0.0);
        W.fill(0.0);
        B.fill(0.0);
        cd.fill(0.0);
        corr.fill(0.0);
        
        for (unsigned int k = 0; k < K; ++k) {
            tmp2 = mu.slice(k)(arma::span(i), arma::span());
            var = var + tmp2.t() * tmp2 
                + qinmatr(sigma.slice(k)(arma::span(i),arma::span())) 
                * weight(k);
            W   = W + qinmatr(sigma.slice(k)(arma::span(i),arma::span()))
                * weight(k);
            d   = arma::trans(mu.slice(k)(arma::span(i), arma::span()))
                - means;
            B   = B + d * d.t() * weight(k);
        }
        var     = var - means * means.t();
        cd      = arma::diagmat(1.0 / arma::sqrt(arma::diagvec(var)));
        corr    = cd * var * cd;
        Rtr     = 1 - arma::trace(W) / arma::trace(var);
        Rdet    = 1 - std::log(arma::det(W)) / std::log(arma::det(var));
        zm.fill(0.0);        
        zm(1)   = 1.0;
        zm(3)   = std::exp(std::log(1.0) + std::log(3.0));        
        for (unsigned int m = 0; m < 4; ++m) {
            for (unsigned int rr = 0; rr < r; ++rr) {
                for (unsigned int k = 0; k < K; ++k) {
                    sigmavec(k) = qinmatr(sigma.slice(k)(arma::span(i), arma::span()))(rr,rr);                    
                }
                tmp = mu.tube(i, rr, i, rr) - means(rr); 
                higher(rr, m)   = arma::as_scalar(weight.t() * arma::pow(tmp, m + 1));
                for (unsigned int n = 0; n < (m + 1); ++n) {
                    arma::vec ss = arma::pow(tmp, m - n);
                    cm = arma::pow(tmp, m - n) % arma::pow(sigmavec, (n + 1) / 2)
                        * zm[n];
                    higher(rr, m)   = higher(rr, m) + R::choose(m + 1, n + 1) 
                        * arma::as_scalar(weight.t() * cm);                    
                }
            }
        }
        skewness    = higher.col(2) / arma::pow(higher.col(1), 1.5);
        kurtosis    = higher.col(3) / arma::pow(higher.col(1), 2);
        meanOut.row(i) = arma::trans(means);
        RtrOut(i)   = Rtr;
        RdetOut(i)  = Rdet;
        for (unsigned int j = 0; j < ij; ++j) {
            corrOut(i, j) = corr(perm(j, 0), perm(j, 1));
        }
        for (unsigned int rr = 0; rr < r; ++rr) {
            varOut(i, rr) = var(rr, rr);
        }
        skewnessOut.row(i) = arma::trans(skewness);
        kurtosisOut.row(i) = arma::trans(kurtosis);
    }

    return Rcpp::List::create(Rcpp::Named("Rtr", RtrOut),
            Rcpp::Named("mean", meanOut),
            Rcpp::Named("Rdet", RdetOut),
            Rcpp::Named("corr", corrOut),
            Rcpp::Named("var", varOut),
            Rcpp::Named("skewness", skewnessOut),
            Rcpp::Named("kurtosis", kurtosisOut));
}

#endif /* __FINMIX_MOMENTS_H__ */



