#include "PriorNormalFix.h"
#include "ParNormalFix.h"
#include "distributions.h"

PriorNormalFix::PriorNormalFix () : HIER(false),
   INDEPENDENT(false)
{
}

PriorNormalFix::PriorNormalFix (const FinmixPrior& prior) :
   HIER(prior.hier),
   INDEPENDENT(prior.type == "condconjugate" ? false : true)
{
   Rcpp::List          tmpMu((SEXP)prior.par["mu"]);
   Rcpp::List          tmpSigma((SEXP)prior.par["sigma"]);
   Rcpp::NumericMatrix tmpb((SEXP)tmpMu["b"]);
   const unsigned int  M = tmpb.nrow();
   const unsigned int  K = tmpb.ncol();

   if (INDEPENDENT)
   {
      Rcpp::NumericMatrix tmpBinv((SEXP)tmpMu["Binv"]);
      BStart = arma::mat(tmpBinv.begin(), M, K, true, true);
      BPost  = BStart;
   }
   else
   {
      Rcpp::NumericMatrix tmpN((SEXP)tmpMu["N0"]);
      BStart = arma::mat(tmpN.begin(), M, K, true, true);
      BPost  = BStart;
   }
   bStart = arma::mat(tmpb.begin(), M, K, true, true);
   Rcpp::NumericMatrix tmpc((SEXP)tmpSigma["c"]);
   Rcpp::NumericMatrix tmpC((SEXP)tmpSigma["C"]);

   cStart = arma::mat(tmpc.begin(), M, K, true, true);
   cPost  = cStart;
   CStart = arma::mat(tmpC.begin(), M, K, true, true);
   CPost  = CStart;
   if (HIER)
   {
      g = tmpSigma["g"];
      G = tmpSigma["G"];
   }
}

inline
void PriorNormalFix::update(const unsigned int& K, const arma::mat& y,
                            arma::ivec& S, const arma::vec& T, ParNormalFix& par)
{
   arma::mat  repY  = arma::repmat(y, 1, K);
   arma::imat repS  = arma::repmat(S, 1, K);
   arma::imat compM = arma::ones<arma::imat>(S.n_elem, K);

   for (unsigned int k = 0; k < K; ++k)
   {
      compM.col(k) = compM.col(k) * (k + 1);
   }
   arma::umat ind       = (repS == compM);
   arma::mat  indDouble = arma::conv_to<arma::mat>::from(ind);

   repY %= indDouble;
   arma::rowvec sprod = sum(repY, 0);
   arma::rowvec sind  = sum(indDouble, 0);

   if (INDEPENDENT)
   {
      if (!par.INDEPENDENT)
      {
         par.INDEPENDENT = true;
      }
      cPost = cStart + 0.5 * sind;
      for (unsigned int k = 0; k < K; ++k)
      {
         CPost(k) = CStart(k);
         arma::uvec yind = find(repY.col(k) != 0.0);
         arma::mat  y    = repY.rows(yind);
         arma::vec  b    = y.col(k) - par.mu(k);
         CPost(k) += 0.5 * arma::as_scalar(arma::trans(b) * b);
      }
      par.sigma = 1.0 / rgammaprod(cPost, CPost);
      arma::rowvec BinvPost = BStart + sind % (1.0 / par.sigma);
      BPost  = 1.0 / BinvPost;
      bPost  = BStart % bStart;
      bPost += 1.0 / par.sigma % sprod;
      bPost %= BPost;
   }
   else      /* conditionally conjugate prior */

   {
      arma::rowvec N0Post = BStart + sind;
      BPost = 1.0 / N0Post;
      bPost = (bStart % BStart + sprod) / N0Post;
      cPost = cStart + 0.5 * sind;
      arma::rowvec ck = BStart % sind / N0Post;
      for (unsigned int k = 0; k < K; ++k)
      {
         if (sind(k) > 0)
         {
            double     yk = sprod(k) / sind(k);
            CPost(k) = CStart(k);
            arma::uvec yind = find(repY.col(k) != 0.0);
            arma::mat  y    = repY.rows(yind);
            arma::vec  sk   = y.col(k) - yk;
            CPost(k) += 0.5 * arma::as_scalar(arma::trans(sk) * sk);
            CPost(k) += 0.5 * (yk - bStart(k)) * (yk - bStart(k)) * ck(k);
         }
         else
         {
            CPost(k)  = CStart(k);
            CPost(k) += 0.5 * (sprod(k) - bStart(k)) * (sprod(k) - bStart(k)) * ck(k);
         }
      }
   }
}

inline
void PriorNormalFix::updateHier(const ParNormalFix& par)
{
   double gN = arma::sum(cStart) + g;
   double GN = arma::sum(1.0 / par.sigma) + G;

   CStart.fill(R::rgamma(gN, 1.0));
   CStart /= GN;
}
