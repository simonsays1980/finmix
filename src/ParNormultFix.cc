#include "ParNormultFix.h"
#include "distributions.h"

ParNormultFix::ParNormultFix (const bool& STARTPAR,
        const FinmixModel& model) : INDEPENDENT(true), 
    mu(model.K, model.r), sigma(model.r, model.r, model.K), 
    sigmainv(model.r, model.r, model.K)
{
    if (model.par.size() > 0) {
        mu = Rcpp::as<arma::mat>(model.par["mu"]);
    }
    if (!STARTPAR && model.K > 1) {
        Rcpp::NumericVector tmpSigma((SEXP) model.par["sigma"]);
        Rcpp::IntegerVector tmpDim = tmpSigma.attr("dim");
        sigma = arma::cube(tmpSigma.begin(), tmpDim[0], tmpDim[1], tmpDim[2], true, true);
    }
}

inline
void ParNormultFix::update (PriorNormultFix& hyperPar) 
{
    if (INDEPENDENT) {
        mu      = rnormult(hyperPar.bPost, hyperPar.BPost);
    } else { /* conditionally conjugate prior */
        for (unsigned int k = 0; k < sigma.n_slices; ++k) {
            sigma.slice(k)      = rinvwishart(hyperPar.cPost(k), 
                    hyperPar.CPost.slice(k));
            sigmainv.slice(k)   = arma::inv(sigma.slice(k));
            hyperPar.BPost.slice(k) = sigma.slice(k) / hyperPar.N0Post(k);            
            hyperPar.BInvPost.slice(k) = arma::inv(hyperPar.BPost.slice(k));
        }       
        mu      = rnormult(hyperPar.bPost, hyperPar.BPost);
    }
}

inline 
void ParNormultFix::permute (const arma::urowvec& compIndex,
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
