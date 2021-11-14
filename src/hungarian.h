/******************************************************************************
*
* Copyright (C) 2013 Lars Simon Zehnder (Bob Pilgrim). All Rights Reserved.
*
* Author: Lars Simon Zehnder <simon.zehnder@gmail.com>
*
* This file is part of the R package finmix.
*
* finmix is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published
* by the Free Software Foundation, either version 3 of the License, or
* any later version.
*
* finmix is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with finmix. If not, see <http://www.gnu.org/licenses/>.
*
******************************************************************************/

#ifndef __FINMIX_HUNGARIAN_H__
#define __FINMIX_HUNGARIAN_H__

#include <RcppArmadillo.h>

/* FORWARD DECLARATION */
void step_one(unsigned int &step, arma::mat &cost,
              const unsigned int &N);

void step_two(unsigned int &step, const arma::mat &cost,
              arma::umat &indM, arma::ivec &rcov,
              arma::ivec &ccov, const unsigned int &N);

void step_three(unsigned int &step, const arma::umat &indM,
                arma::ivec &ccov, const unsigned int &N);

void step_four(unsigned int &step, const arma::mat &cost,
               arma::umat &indM, arma::ivec &rcov,
               arma::ivec &ccov, int &rpath_0, int &cpath_0,
               const unsigned int &N);

void step_five(unsigned int &step, arma::umat &indM,
               arma::ivec &rcov, arma::ivec &ccov,
               arma::imat &path, int &rpath_0, int &cpath_0,
               const unsigned int &N);

void step_six(unsigned int &step, arma::mat &cost,
              const arma::ivec &rcov, const arma::ivec &ccov,
              const unsigned int &N);

void find_noncovered_zero(int &row, int &col,
                          const arma::mat &cost, const arma::ivec &rcov,
                          const arma::ivec &ccov, const unsigned int &N);

bool star_in_row(int &row, const arma::umat &indM,
                 const unsigned int &N);

void find_star_in_col(const int &col, int &row,
                      const arma::umat &indM, const unsigned int &N);

void find_star_in_row(const int &row, int &col,
                      const arma::umat &indM, const unsigned int &N);

void find_prime_in_row(const int &row, int &col,
                       const arma::umat &indM, const unsigned int &N);

void augment_path(const unsigned int &path_count, arma::umat &indM,
                  const arma::imat &path);

void clear_covers(arma::ivec &rcov, arma::ivec &ccov);

void erase_primes(arma::umat &indM, const unsigned int &N);

void find_smallest(double &minval, const arma::mat &cost,
                   const arma::ivec &rcov, const arma::ivec &ccov,
                   const unsigned int &N);

/* ALGORITHM */

/**
 * Searches for the assignment of rows to columns
 * with the minimum cost regarding the cost matrix 'cost'.
 * To find a maximal assignment the cost matrix can be
 * transformed via DBL_MAX - cost.
 * @param input_cost    constant Armadillo matrix reference.
 * @see step_one(), step_two(), step_three(), step_four(),
 *      step_five(), step_six()
 * @return uword Armadillo matrix.
 * */
inline
arma::umat hungarian(const arma::mat &input_cost)
{
   const unsigned int N       = input_cost.n_rows;
   unsigned int       step    = 1;
   int                cpath_0 = 0;
   int                rpath_0 = 0;
   arma::mat          cost(input_cost);
   arma::umat         indM(N, N);
   arma::ivec         rcov(N);
   arma::ivec         ccov(N);
   arma::imat         path(2 * N, 2);

   indM = arma::zeros<arma::umat>(N, N);
   bool done = false;

   while (!done)
   {
      switch (step)
      {
      case 1:
         step_one(step, cost, N);
         break;

      case 2:
         step_two(step, cost, indM, rcov, ccov, N);
         break;

      case 3:
         step_three(step, indM, ccov, N);
         break;

      case 4:
         step_four(step, cost, indM, rcov, ccov,
                   rpath_0, cpath_0, N);
         break;

      case 5:
         step_five(step, indM, rcov, ccov,
                   path, rpath_0, cpath_0, N);
         break;

      case 6:
         step_six(step, cost, rcov, ccov, N);
         break;

      case 7:
         done = true;
         break;
      }
   }
   return indM;
}

/**
 * Reduce each row by its minimum.
 * @param step  unsigned integer reference.
 * @param cost  Armadillo matrix reference.
 * @param N     constant unsigned integer reference.
 * @see step_two(), hungarian()
 * @return void
 * */
inline
void step_one(unsigned int &step, arma::mat &cost,
              const unsigned int &N)
{
   for (unsigned int r = 0; r < N; ++r)
   {
      cost.row(r) -= arma::min(cost.row(r));
   }
   step = 2;
}

/**
 * Find a zero in the resulting cost matrix of step one.
 * @param step  unsigned integer reference.
 * @param cost  constant Armadillo matrix reference.
 * @param indM  uword Armadillo matrix reference.
 * @param rcov  iword Armadillo vector reference.
 * @param ccov  iword Armadillo vector reference.
 * @param N     constant unsigned integer reference.
 * @see step_three(), step_one()
 * @return void
 * */
inline
void step_two(unsigned int &step, const arma::mat &cost,
              arma::umat &indM, arma::ivec &rcov,
              arma::ivec &ccov, const unsigned int &N)
{
   for (unsigned int r = 0; r < N; ++r)
   {
      for (unsigned int c = 0; c < N; ++c)
      {
         if (cost.at(r, c) == 0.0 && rcov.at(r) == 0 && ccov.at(c) == 0)
         {
            indM.at(r, c) = 1;
            rcov.at(r)    = 1;
            ccov.at(c)    = 1;
            break;                                                  // Only take the first
                                                                    // zero in a row and column
         }
      }
   }
   /* for later reuse */
   rcov.fill(0);
   ccov.fill(0);
   step = 3;
}

/**
 * Cover each column containing a starred zero. If N
 * columns are covered the starred zeros describe a
 * complete set of unqiue assignments. In this case
 * go to the last STEP 7 (hungarian()) otherwise go
 * to STEP 4 (step_four()).
 * @param step  unsigned integer reference.
 * @param indM  constant uword Armadillo matrix reference.
 * @param rcov  iword Armadillo vector reference.
 * @param ccov  iword Armadillo vector reference.
 * @param N     constant unsigned integer reference.
 * @see hungarian(), step_four()
 * @return void
 * */
inline
void step_three(unsigned int &step, const arma::umat &indM,
                arma::ivec &ccov, const unsigned int &N)
{
   unsigned int colcount = 0;

   for (unsigned int r = 0; r < N; ++r)
   {
      for (unsigned int c = 0; c < N; ++c)
      {
         if (indM.at(r, c) == 1)
         {
            ccov.at(c) = 1;
         }
      }
   }
   for (unsigned int c = 0; c < N; ++c)
   {
      if (ccov.at(c) == 1)
      {
         ++colcount;
      }
   }
   if (colcount == N)
   {
      step = 7;
   }
   else
   {
      step = 4;
   }
}

/**
 * Find a noncovered zero and prime it. If there
 * is no starred zero in the row containing this
 * primed zero. Go to STEP 5 (step_five()). Otherwise,
 * cover this row and uncover the column containing the
 * starred zero. Continue this way until there are
 * no uncovered zeros left. Save the smallest
 * uncovered VALUE (not necessary a zero) and go
 * to STEP 6 (step_six()).
 * @param step      unsigned integer reference.
 * @param cost      const Armadillo matrix reference.
 * @param indM      uword Armadillo matrix reference.
 * @param rcov      iword Armadillo vector reference.
 * @param ccov      iword Armadillo vector reference.
 * @param rpath_0   integer reference.
 * @param cpath_0   integer reference.
 * @param N         constant unsigned integer reference.
 * @see step_five(), step_six(), star_in_row(),
 *      find_star_in_row()
 * @return void
 **/
inline
void step_four(unsigned int &step, const arma::mat &cost,
               arma::umat &indM, arma::ivec &rcov, arma::ivec &ccov,
               int &rpath_0, int &cpath_0, const unsigned int &N)
{
   int  row  = -1;
   int  col  = -1;
   bool done = false;

   while (!done)
   {
      find_noncovered_zero(row, col, cost, rcov,
                           ccov, N);

      if (row == -1)
      {
         done = true;
         step = 6;
      }
      else
      {
         /* uncovered zero */
         indM(row, col) = 2;
         if (star_in_row(row, indM, N))
         {
            find_star_in_row(row, col, indM, N);
            /* Cover the row with the starred zero
             * and uncover the column with the starred
             * zero.
             */
            rcov.at(row) = 1;
            ccov.at(col) = 0;
         }
         else
         {
            /* No starred zero in row with
             * uncovered zero
             */
            done    = true;
            step    = 5;
            rpath_0 = row;
            cpath_0 = col;
         }
      }
   }
}

/**
 * Construct a series of alternating primed and starred
 * zeros as follows. Let Z0 represent the uncovered primed
 * zero found in Step 4. Let Z1 denote the starred zero in
 * the column of Z0 (if any). Let Z2 denote the primed zero
 * in the row of Z1 (there will always be one, given that
 * there is a Z1). Continue until the series terminates
 * at a primed zero that has no starred zero in its
 * column. Unstar each starred zero of the series, star
 * each primed zero of the series, erase all primes and
 * uncover every line in the matrix. Return to STEP 3
 * (step_three()).
 * @param step      unsigned integer reference.
 * @param indM      uword Armadillo matrix reference.
 * @param rcov      iword Armadillo vector reference.
 * @param ccov      iword Armadillo vector reference.
 * @param path      iword Armadillo matrix reference.
 * @param rpath_0   integer reference.
 * @param cpath_0   integer reference.
 * @param N         constant unsigned integer reference.
 * @see step_four(), step_three()
 * @return void
 * */
inline
void step_five(unsigned int &step,
               arma::umat &indM, arma::ivec &rcov,
               arma::ivec &ccov, arma::imat &path,
               int &rpath_0, int &cpath_0,
               const unsigned int &N)
{
   bool         done       = false;
   int          row        = -1;
   int          col        = -1;
   unsigned int path_count = 1;

   path.at(path_count - 1, 0) = rpath_0;
   path.at(path_count - 1, 1) = cpath_0;
   while (!done)
   {
      find_star_in_col(path.at(path_count - 1, 1), row,
                       indM, N);
      if (row > -1)
      {
         /* Starred zero in row 'row' */
         ++path_count;
         path.at(path_count - 1, 0) = row;
         path.at(path_count - 1, 1) = path.at(path_count - 2, 1);
      }
      else
      {
         done = true;
      }
      if (!done)
      {
         /* If there is a starred zero find a primed
         * zero in this row; write index to 'col' */
         find_prime_in_row(path.at(path_count - 1, 0), col,
                           indM, N);
         ++path_count;
         path.at(path_count - 1, 0) = path.at(path_count - 2, 0);
         path.at(path_count - 1, 1) = col;
      }
   }
   augment_path(path_count, indM, path);
   clear_covers(rcov, ccov);
   erase_primes(indM, N);
   step = 3;
}

/**
 * Adds the VALUE (not necessary zero) found in STEP 4
 * (step_four()) to every element of each covered row,
 * and subtracts it from every element of each uncovered
 * column. Returns to STEP 4 (step_four()) without
 * altering any starred or primed zeros nor covered lines.
 * @param step  unsigned integer reference.
 * @param cost  Armadillo matrix reference.
 * @param rcov  constant iword Armadillo vector reference.
 * @param ccov  constant iword Armadillo vector reference.
 * @param N     constant unsigned integer reference.
 * @see step_four(), find_smallest()
 * @return void
 * */
inline
void step_six(unsigned int &step, arma::mat &cost,
              const arma::ivec &rcov, const arma::ivec &ccov,
              const unsigned int &N)
{
   double minval = DBL_MAX;

   find_smallest(minval, cost, rcov, ccov, N);
   for (unsigned int r = 0; r < N; ++r)
   {
      for (unsigned int c = 0; c < N; ++c)
      {
         if (rcov.at(r) == 1)
         {
            cost.at(r, c) += minval;
         }
         if (ccov.at(c) == 0)
         {
            cost.at(r, c) -= minval;
         }
      }
   }
   step = 4;
}

/* Helper functions */

/**
 * Finds noncovered zeros in the cost matrix.
 * @param row   integer reference.
 * @param col   integer reference.
 * @param cost  constant Armadillo matrix reference.
 * @param rcov  constant iword Armadillo vector reference.
 * @param ccov  constant iword Armadillo vector reference.
 * @param N     constant unsigned integer reference.
 * @see step_four()
 * @return void
 */
inline
void find_noncovered_zero(int &row, int &col,
                          const arma::mat &cost, const arma::ivec &rcov,
                          const arma::ivec &ccov, const unsigned int &N)
{
   unsigned int r = 0;
   unsigned int c;
   bool         done = false;

   row = -1;
   col = -1;
   while (!done)
   {
      c = 0;
      while (true)
      {
         if (cost.at(r, c) == 0.0 && rcov.at(r) == 0 && ccov.at(c) == 0)
         {
            row  = r;
            col  = c;
            done = true;
         }
         ++c;
         if (c == N || done)
         {
            break;
         }
      }
      ++r;
      if (r == N)
      {
         done = true;
      }
   }
}

/**
 * Indicates if a starred zero is contained in a certain
 * row of the cost matrix searching the indicator matrix.
 * @param row   integer reference.
 * @param indM  constant uword Armadillo matrix reference.
 * @param N     constant unsigned integer reference.
 * @see step_four()
 * @return A bool indicating if there is.
 * */
inline
bool star_in_row(int &row, const arma::umat &indM,
                 const unsigned int &N)
{
   bool tmp = false;

   for (unsigned int c = 0; c < N; ++c)
   {
      if (indM.at(row, c) == 1)
      {
         tmp = true;
         break;
      }
   }
   return tmp;
}

/**
 * Finds a starred zero in a certain column of the cost
 * matrix by searching the indicator matrix and writing
 * to the references 'row'.
 * @param col   constant integer reference.
 * @param row   integer reference.
 * @param indM  constant uword Armadillo matrix reference.
 * @param N     constant unsigned integer reference.
 * @see step_five()
 * @return void
 * */
inline
void find_star_in_col(const int &col, int &row,
                      const arma::umat &indM, const unsigned int &N)
{
   row = -1;
   for (unsigned int r = 0; r < N; ++r)
   {
      if (indM.at(r, col) == 1)
      {
         row = r;
      }
   }
}

/**
 * Finds a starred zero in a certain row of the cost
 * matrix by searching the indicator matrix and writing
 * to the references 'col'.
 * @param row   constant integer reference.
 * @param col   integer reference.
 * @param indM  constant uword Armadillo matrixi reference.
 * @param N     constant unsigned integer reference.
 * @see step_four()
 * @return void
 * */
inline
void find_star_in_row(const int &row, int &col,
                      const arma::umat &indM, const unsigned int &N)
{
   col = -1;
   for (unsigned int c = 0; c < N; ++c)
   {
      if (indM.at(row, c) == 1)
      {
         col = c;
      }
   }
}

/**
 * Finds a primed zero in a certain row of the cost
 * matrix by searching the indicator matrix. Writes
 * result to the argument 'col'.
 * @param row   constant integer reference.
 * @param col   integer reference.
 * @param indM  constant uword Armadillo matrix reference.
 * @param N     constant unsigned integer reference.
 * @see step_five()
 * @return void
 * */
inline
void find_prime_in_row(const int &row, int &col,
                       const arma::umat &indM, const unsigned int &N)
{
   for (unsigned int c = 0; c < N; ++c)
   {
      if (indM.at(row, c) == 2)
      {
         col = c;
      }
   }
}

/**
 * Augments the path through the cost matrix starting at
 * a primed zero (found in step_four()) and ends at a
 * primed zero (see bipartite graph theory).
 * @param path_count    constant integer reference.
 * @param indM          uword Armadillo matrix reference.
 * @param path          constant iword Armadillo matrix reference.
 * @see step_five()
 * @return void
 * */
inline
void augment_path(const unsigned int &path_count, arma::umat &indM,
                  const arma::imat &path)
{
   for (unsigned int p = 0; p < path_count; ++p)
   {
      if (indM.at(path(p, 0), path(p, 1)) == 1)
      {
         indM.at(path(p, 0), path(p, 1)) = 0;
      }
      else
      {
         indM.at(path(p, 0), path(p, 1)) = 1;
      }
   }
}

/**
 * Clears the cover vectors.
 * @param rcov  iword Armadillo vector reference.
 * @param ccov  iword Armadillo vector reference.
 * @see step_five()
 * @return void
 * */
inline
void clear_covers(arma::ivec &rcov, arma::ivec &ccov)
{
   rcov.fill(0);
   ccov.fill(0);
}

/**
 * Erases all indicated primed zeros from the indicator
 * matrix.
 * @param indM  uword Armadillo matrix reference.
 * @param N     constant unsigned integer reference.
 * @see step_five()
 * @return void.
 * */
inline
void erase_primes(arma::umat &indM, const unsigned int &N)
{
   for (unsigned int r = 0; r < N; ++r)
   {
      for (unsigned int c = 0; c < N; ++c)
      {
         if (indM.at(r, c) == 2)
         {
            indM.at(r, c) = 0;
         }
      }
   }
}

/**
 * Finds smallest value in the cost matrix over all
 * uncovered columns and rows and writes it to
 * 'minval'
 * @param minval    double reference.
 * @param cost      constant Armadillo matrix reference.
 * @param rcov      constant iword Amradillo vector reference.
 * @param ccov      constant iword Armadillo vector reference.
 * @param N         constant unsigned integer reference.
 * @see step_six()
 * @return void
 * */
inline
void find_smallest(double &minval, const arma::mat &cost,
                   const arma::ivec &rcov, const arma::ivec &ccov,
                   const unsigned int &N)
{
   for (unsigned int r = 0; r < N; ++r)
   {
      for (unsigned int c = 0; c < N; ++c)
      {
         if (rcov.at(r) == 0 && ccov.at(c) == 0)
         {
            if (minval > cost.at(r, c))
            {
               minval = cost.at(r, c);
            }
         }
      }
   }
}

#endif /* __FINMIX_HUNGARIAN_H__ */
