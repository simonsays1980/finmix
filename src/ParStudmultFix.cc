#include "ParStudmultFix.h"
#include "distributions.h"

ParStudmultFix::ParStudmultFix (const bool& STARTPAR,
        const FinmixModel& model) : mu(model.K, model.r), 
        sigma(model.r, model.r, model.K), 
        sigmainv(model.r, model.r, model.K), df(model.K), 
        acc(model.K), INDEPENDENT(true)
{
    acc.fill(0.0);
    if (model.par.size() > 0) {
        if (!Rf_isNull(model.par["mu"])) {
            mu = Rcpp::as<arma::mat>(model.par["mu"]);
        }
        if (!Rf_isNull(model.par["df"])) {
            df = Rcpp::as<arma::rowvec>(model.par["df"]);
        }
    }
    if (!STARTPAR && model.K > 1) {
        Rcpp::NumericVector tmpSigma((SEXP) model.par["sigma"]);
        Rcpp::IntegerVector tmpDim = tmpSigma.attr("dim");
        sigma = arma::cube(tmpSigma.begin(), tmpDim[0], tmpDim[1], tmpDim[2], true, true);
        Rcpp::NumericVector tmpDf((SEXP) model.par["df"]);
        df = arma::rowvec(tmpDf.begin(), model.K, true, true);
    }
}

inline
void ParStudmultFix::update (PriorStudmultFix& hyperPar) 
{
    /* See PriorStudmultFix.cc */
}

inline 
void ParStudmultFix::permute (const arma::urowvec& compIndex,
        const arma::urowvec& permIndex) 
{
    mu.cols(compIndex) = mu.cols(permIndex);
    arma::cube tmpSigma = sigma;
    arma::cube tmpSigmaInv = sigmainv;
    for (unsigned int k = 0; k < sigma.n_slices; ++k) {
        sigma.slice(compIndex(k)) = tmpSigma.slice(permIndex(k));
        sigmainv.slice(compIndex(k)) = tmpSigmaInv.slice(permIndex(k));
    }    
}
