#include "ParStudentFix.h"

ParStudentFix::ParStudentFix (const bool& STARTPAR,
                              const FinmixModel& model) : mu(model.K),
   sigma(model.K), acc(model.K), INDEPENDENT(false)
{
   if (!Rf_isNull(model.par))
   {
      if (!Rf_isNull(model.par["mu"]))
      {
         mu = Rcpp::as<arma::rowvec>(model.par["mu"]);
      }
      if (!Rf_isNull(model.par["df"]))
      {
         df = Rcpp::as<arma::rowvec>(model.par["df"]);
      }
   }
   if (!STARTPAR && model.K > 1)
   {
      arma::rowvec tmpsigma = Rcpp::as<arma::rowvec>
                                 ((SEXP)model.par["sigma"]);
      sigma = tmpsigma;
   }
}

inline
void ParStudentFix::update(const PriorStudentFix& hyperPar)
{
   /* See PriorStudentFix.cc */
}

inline
void ParStudentFix::permute(const arma::urowvec& compIndex,
                            const arma::urowvec& permIndex)
{
   mu(compIndex)    = mu(permIndex);
   sigma(compIndex) = sigma(permIndex);
   df(compIndex)    = df(permIndex);
}
