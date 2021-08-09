#include "ParStudmultInd.h"
#include "distributions.h"

ParStudmultInd::ParStudmultInd (const bool& STARTPAR,
        const FinmixModel& model) :
    ParStudmultFix(STARTPAR, model),
    weight(model.K)
{
    if (!STARTPAR && model.K > 1) {
        weight = model.weight;
    }
}

void ParStudmultInd::update (PriorStudmultInd& hyperPar) 
{
    ParStudmultFix::update(hyperPar);
    weight = rdirichlet(hyperPar.weightPost);
}
