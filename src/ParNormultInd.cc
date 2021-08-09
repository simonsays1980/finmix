#include "ParNormultInd.h"
#include "distributions.h"

ParNormultInd::ParNormultInd (const bool& STARTPAR,
        const FinmixModel& model) :
    ParNormultFix(STARTPAR, model),
    weight(model.K)
{
    if (!STARTPAR && model.K > 1) {
        weight = model.weight;
    }
}

void ParNormultInd::update (PriorNormultInd& hyperPar) 
{
    ParNormultFix::update(hyperPar);
    weight = rdirichlet(hyperPar.weightPost);
}
