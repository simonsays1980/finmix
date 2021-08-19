#include "ParNormalFix.h"
#include "distributions.h"

ParNormalFix::ParNormalFix (const bool& STARTPAR,
        const FinmixModel& model) : mu(model.K), sigma(model.K),
    INDEPENDENT(false)
{ 
    if (!Rf_isNull(model.par)) {
        mu = Rcpp::as<arma::rowvec>(model.par["mu"]);
    }
    if (!STARTPAR && model.K > 1) {
        sigma = Rcpp::as<arma::rowvec>(model.par["sigma"]);
    }
}

void ParNormalFix::update (const PriorNormalFix& hyperPar)
{
    if (INDEPENDENT) {
        mu      = rnormal(hyperPar.bPost, hyperPar.BPost);       
    } else { /* conditionally conjugate prior */
        sigma   = 1.0 / rgammaprod(hyperPar.cPost, hyperPar.CPost); 
        mu      = rnormal(hyperPar.bPost, sigma % hyperPar.BPost);
    }
}

void ParNormalFix::permute (const arma::urowvec& compIndex,
        const arma::urowvec& permIndex)
{
    mu(compIndex)       = mu(permIndex);
    sigma(compIndex)    = sigma(permIndex);
}


