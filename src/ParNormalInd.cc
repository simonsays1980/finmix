#include "ParNormalInd.h"
#include "distributions.h"

ParNormalInd::ParNormalInd (const bool& STARTPAR,
        const FinmixModel& model) :
    ParNormalFix(STARTPAR, model),
    weight(model.K)
{
    if (!STARTPAR && model.K > 1) {
        weight = model.weight;
    }
}

void ParNormalInd::update (const PriorNormalInd& hyperPar)
{
    ParNormalFix::update(hyperPar);
    weight = rdirichlet(hyperPar.weightPost);
}
