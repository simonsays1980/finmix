#include "PriorStudmultInd.h"
#include "ParStudmultInd.h"
#include "posterior.h"

PriorStudmultInd::PriorStudmultInd (const FinmixPrior& prior) :
   PriorStudmultFix(prior), weightStart(prior.weight),
   weightPost(prior.weight)
{
}

inline
void PriorStudmultInd::update(const unsigned int& K, const arma::mat& y,
                              arma::ivec& S, const arma::vec& T, ParStudmultInd& par)
{
   PriorStudmultFix::update(K, y, S, T, par);
   weightPost = posterior_multinomial(K, S, weightStart);
}
