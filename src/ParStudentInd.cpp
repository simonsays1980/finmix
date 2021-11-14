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

void ParStudentInd::update(const PriorStudentInd& hyperPar)
{
   // updating the parameters is performed in PriorStudentInd.update()
   weight = rdirichlet(hyperPar.weightPost);
}
