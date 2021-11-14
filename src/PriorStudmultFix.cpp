#include "PriorStudmultFix.h"
#include "ParStudmultFix.h"
#include "distributions.h"
#include "DataClass.h"
#include "likelihood.h"
#include "prior_likelihood.h"

PriorStudmultFix::PriorStudmultFix () : HIER(false),
   INDEPENDENT(false)
{
}

PriorStudmultFix::PriorStudmultFix (const FinmixPrior& prior) :
   HIER(prior.hier),
   INDEPENDENT(prior.type == "condconjugate" ? false : true)
{
   Rcpp::List          tmpMu((SEXP)prior.par["mu"]);
   Rcpp::List          tmpSigma((SEXP)prior.par["sigma"]);
   Rcpp::NumericMatrix tmpb((SEXP)tmpMu["b"]);
   const unsigned int  r = tmpb.nrow();
   const unsigned int  K = tmpb.ncol();

   if (INDEPENDENT)
   {
      Rcpp::NumericVector tmpBinv((SEXP)tmpMu["Binv"]);
      Rcpp::IntegerVector tmpDims = tmpBinv.attr("dim");
      BInvStart = arma::cube(tmpBinv.begin(), tmpDims[0], tmpDims[1],
                             tmpDims[2], true, true);
      BStart = arma::cube(tmpDims[0], tmpDims[1], tmpDims[2]);
      for (unsigned int k = 0; k < K; ++k)
      {
         BStart.slice(k) = arma::inv(BInvStart.slice(k));
      }
      BPost    = BStart;
      BInvPost = BInvStart;
      /* Initialize N0Start and N0Post (not used if INDEPENDENT) */
      N0Start = arma::rowvec(r);
      N0Post  = N0Start;
   }
   else      /* conditionally conjugate prior */
   {
      Rcpp::NumericVector tmpN((SEXP)tmpMu["N0"]);
      N0Start   = arma::rowvec(tmpN.begin(), K, true, true);
      N0Post    = N0Start;
      BStart    = arma::cube(r, r, K);
      BInvStart = BStart;
      BPost     = BStart;
      BInvPost  = BInvStart;
   }
   bStart = arma::mat(tmpb.begin(), r, K, true, true);
   bPost  = bStart;
   Rcpp::NumericMatrix tmpc((SEXP)tmpSigma["c"]);
   Rcpp::NumericVector tmpC((SEXP)tmpSigma["C"]);
   Rcpp::IntegerVector tmpDims2 = tmpC.attr("dim");

   cStart = arma::rowvec(tmpc.begin(), K, true, true);
   cPost  = cStart;
   CStart = arma::cube(tmpC.begin(), tmpDims2[0],
                       tmpDims2[1], tmpDims2[2], true, true);
   CPost = CStart;
   Rcpp::NumericVector tmpLogdetC((SEXP)tmpSigma["logdetC"]);

   logdetC = arma::rowvec(tmpLogdetC.begin(), K, true, true);
   if (HIER)
   {
      g = tmpSigma["g"];
      Rcpp::NumericMatrix tmpG((SEXP)tmpSigma["G"]);
      G = arma::mat(tmpG.begin(), r, r, true, true);
   }
   Rcpp::List tmpDf((SEXP)prior.par["df"]);

   dftype = Rcpp::as<std::string>(tmpDf["type"]);
   trans  = Rcpp::as<double>(tmpDf["trans"]);
   a0     = Rcpp::as<double>(tmpDf["a0"]);
   b0     = Rcpp::as<double>(tmpDf["b0"]);
   d      = Rcpp::as<double>(tmpDf["d"]);
   Rcpp::NumericVector tmpMhTune((SEXP)tmpDf["mhtune"]);

   mhTune = arma::rowvec(tmpMhTune.begin(), K, true, true);
}

void PriorStudmultFix::update(const unsigned int& K, const arma::mat& y,
                              arma::ivec& S, const arma::vec& T, ParStudmultFix& par)
{
   arma::mat  repY  = arma::repmat(y, 1, K);
   arma::imat repS  = arma::repmat(S, 1, K);
   arma::imat compM = arma::ones<arma::imat>(S.n_elem, K);

   for (unsigned int k = 0; k < K; ++k)
   {
      compM.col(k) = compM.col(k) * (k + 1);
   }
   arma::umat   ind       = (repS == compM);
   arma::mat    indDouble = arma::conv_to<arma::mat>::from(ind);
   arma::rowvec sind      = sum(indDouble, 0);

   if (INDEPENDENT)
   {
      if (!par.INDEPENDENT)
      {
         par.INDEPENDENT = true;
      }
      cPost = cStart + 0.5 * sind;
      double sign = 0.0;
      for (unsigned int k = 0; k < K; ++k)
      {
         CPost.slice(k) = CStart.slice(k);
         arma::uvec yind = find(ind.col(k) != 0.0);
         arma::mat  y2   = y.rows(yind);
         arma::mat  b    = y2;
         b.each_row()   -= arma::trans(par.mu.col(k));
         CPost.slice(k) += 0.5 * arma::trans(b) * b;
         arma::log_det(logdetC(k), sign, CPost.slice(k));
         logdetC(k)         = logdetC(k) * sign;
         par.sigma.slice(k) = rinvwishart(cPost(k), CPost.slice(k));
         par.sigma.slice(k);
         par.sigmainv.slice(k) = arma::inv(par.sigma.slice(k));
         BInvPost.slice(k)     = BInvStart.slice(k) + sind(k) * arma::inv(par.sigma.slice(k));
         BPost.slice(k)        = arma::inv(BInvPost.slice(k));
         bPost.col(k)          = BInvStart.slice(k) * bStart.col(k)
                                 + arma::inv(par.sigma.slice(k)) * arma::trans(arma::sum(y2, 0));
         bPost.col(k) = BPost.slice(k) * bPost.col(k);
      }
   }
   else      /* conditionally conjugate prior */
   {
      if (par.INDEPENDENT)
      {
         par.INDEPENDENT = false;
      }
      /* BStart is actually N0Start */
      N0Post = N0Start + sind;
      for (unsigned int k = 0; k < K; ++k)
      {
         arma::uvec yind = find(ind.col(k) != 0.0);
         arma::mat  y2   = y.rows(yind);
         bPost.col(k)  = bStart.col(k) * N0Start(k);
         bPost.col(k) += arma::trans(arma::sum(y2, 0));
         bPost.col(k) /= N0Post(k);
      }
      double       sign = 0.0;
      cPost = cStart + 0.5 * sind;
      arma::rowvec ck = N0Start % sind / N0Post;
      for (unsigned int k = 0; k < K; ++k)
      {
         CPost.slice(k) = CStart.slice(k);
         arma::uvec   yind = find(ind.col(k) != 0.0);
         arma::mat    y2   = y.rows(yind);
         arma::rowvec sk   = arma::sum(y2, 0);
         if (sind(k) > 0)
         {
            arma::rowvec yk = sk / sind(k);
            arma::mat    dk = y2;
            dk.each_row()  -= yk;
            CPost.slice(k) += 0.5 * arma::trans(dk) * dk;
            CPost.slice(k) += 0.5 * arma::trans(yk - arma::trans(bStart.col(k)))
                              * (yk - arma::trans(bStart.col(k))) * ck(k);
         }
         else
         {
            CPost.slice(k) += 0.5 * (arma::trans(sk) - bStart.col(k))
                              * (arma::trans(sk) - bStart.col(k)) * ck(k);
         }
         arma::log_det(logdetC(k), sign, CPost.slice(k));
         logdetC(k) = logdetC(k) * sign;
      }
   }

   /* The parameter update is done here, as we need data 'y'
    * and classifications 'S' for the Metropolis-Hastings
    * algorithm for the degrees of freedoms update */
   if (par.INDEPENDENT)
   {
      par.mu = rnormult(bPost, BPost);
   }
   else      /* conditionally conjugate prior */
   {
      for (unsigned int k = 0; k < par.sigma.n_slices; ++k)
      {
         par.sigma.slice(k)    = rinvwishart(cPost(k), CPost.slice(k));
         par.sigmainv.slice(k) = arma::inv(par.sigma.slice(k));
         BPost.slice(k)        = par.sigma.slice(k) / N0Post(k);
         BInvPost.slice(k)     = arma::inv(BPost.slice(k));
      }
      par.mu = rnormult(bPost, BPost);
   }
   updateDf(K, y, S, par);
}

void PriorStudmultFix::updateDf(const unsigned int& K,
                                const arma::mat& y, const arma::ivec& S,
                                ParStudmultFix& par)
{
   arma::rowvec dfnew = par.df;
   liklist      lik   = likelihood_studmult(y, par.mu, par.sigmainv,
                                            par.df);
   DataClass    dataC;
   double       loglik   = 0.0;
   double       priorlik = priormixlik_studmult(INDEPENDENT, HIER,
                                                bStart, BInvStart, BStart, cStart, CStart, logdetC, g, G,
                                                par.mu, par.sigma, par.df, trans, a0, b0, d);

   dataC  = classification_fix(K, S, lik);
   loglik = arma::sum(dataC.logLikCd);
   double         urnd        = 0.0;
   double         logliknew   = 0.0;
   double         priorliknew = 0.0;
   double         acc         = 0.0;
   Rcpp::RNGScope scope;

   for (unsigned int k = 0; k < K; ++k)
   {
      urnd     = mhTune(k) * (2.0 * R::runif(0.0, 1.0) - 1.0);
      dfnew(k) = trans + (par.df(k) - trans) * std::exp(urnd);
      liklist lik2 = likelihood_studmult(y, par.mu, par.sigmainv, dfnew);
      dataC       = classification_fix(K, S, lik2);
      logliknew   = arma::sum(dataC.logLikCd);
      priorliknew = priormixlik_studmult(INDEPENDENT, HIER,
                                         bStart, BInvStart, BStart, cStart, CStart, logdetC, g, G,
                                         par.mu, par.sigma, dfnew, trans, a0, b0, d);
      acc = logliknew + priorliknew - (loglik + priorlik) + urnd;
      if (std::log(R::runif(0.0, 1.0)) < acc)
      {
         par.df(k)  = dfnew(k);
         loglik     = logliknew;
         priorlik   = priorliknew;
         par.acc(k) = 1.0;
      }
   }
}

void PriorStudmultFix::updateHier(const ParStudmultFix& par)
{
   if (HIER)
   {
      const unsigned int K  = cPost.n_elem;
      double             gN = arma::sum(cStart) + g;
      arma::mat          GN = G;
      for (unsigned int k = 0; k < K; ++k)
      {
         GN += par.sigma.slice(k);
      }
      CStart.slice(0) = arma::inv(rinvwishart(gN, GN));
      if (K > 1)
      {
         for (unsigned int k = 1; k < K; ++k)
         {
            CStart.slice(k) = CStart.slice(0);
         }
      }
   }
}
