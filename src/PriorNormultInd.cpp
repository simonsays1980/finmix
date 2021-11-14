#include "PriorNormultInd.h"
#include "ParNormultInd.h"
#include "posterior.h"

PriorNormultInd::PriorNormultInd (const FinmixPrior& prior) :
   PriorNormultFix(prior), weightStart(prior.weight),
   weightPost(prior.weight)
{
}

void PriorNormultInd::update(const unsigned int& K, const arma::mat& y,
                             arma::ivec& S, const arma::vec& T, ParNormultInd& par)
{
   PriorNormultFix::update(K, y, S, T, par);
   weightPost = posterior_multinomial(K, S, weightStart);
}
