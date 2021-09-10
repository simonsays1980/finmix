#include "ParStudentInd.h"
#include "distributions.h"

ParStudentInd::ParStudentInd (const bool& STARTPAR,
                              const FinmixModel& model) :
   ParStudentFix(STARTPAR, model),
   weight(model.K)
{
   if (!STARTPAR && model.K > 1)
   {
      weight = model.weight;
   }
}

inline
void ParStudentInd::update(const PriorStudentInd& hyperPar)
{
   weight = rdirichlet(hyperPar.weightPost);
}
