#include "PriorNormalInd.h"
#include "ParNormalInd.h"
#include "posterior.h"

PriorNormalInd::PriorNormalInd (const FinmixPrior& prior) :
    PriorNormalFix(prior),
    weightStart(prior.weight),
    weightPost(prior.weight) {}

inline
void PriorNormalInd::update (const unsigned int& K, const arma::mat& y,
        arma::ivec& S, const arma::vec& T, ParNormalInd& par)
{
    PriorNormalFix::update(K, y, S, T, par);
    weightPost = posterior_multinomial(K, S, weightStart);
}   
