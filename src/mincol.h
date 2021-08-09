/******************************************************************************
 *
 * TODO: Project Title
 *
 * Copyright (C) 2003-2009 ascolab GmbH. All Rights Reserved.
 * Web: http://www.ascolab.com
 *
 * Author: Gerhard Gappmeier <gerhard.gappmeier@ascolab.com>
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 ******************************************************************************/

#ifndef __FINMIX_MINCOL_H__
#define __FINMIX_MINCOL_H__

inline
arma::vec mincol (const arma::mat& m) 
{
    const unsigned int r = m.n_rows;
    arma::vec output(r * (r + 1) / 2);
    unsigned int index = 0;
    for (unsigned int rr = 0; rr < r; ++rr) {
        output.rows(index, index + rr) = m(arma::span(0, rr), rr);
        index = index + rr + 1;
    }
    return output;
}

inline
arma::rowvec minrow (const arma::mat& m)
{
    const unsigned int r = m.n_rows;
    arma::rowvec output(r * (r + 1) / 2);
    unsigned int index = 0;
    for (unsigned int rr = 0; rr < r; ++rr) {
        output.cols(index, index + rr) = arma::trans(m(arma::span(0, rr), rr));
        index = index + rr + 1;
    }
    return output;
}

inline
arma::mat cincolmat (const arma::cube& c) 
{
    const unsigned int r = c.n_rows;
    const unsigned int K = c.n_slices;
    arma::mat output(r * (r + 1) / 2, K);
    unsigned int index = 0;
    for (unsigned int k = 0; k < K; ++k) {
        for (unsigned int rr = 0; rr < r; ++rr) {
            output(arma::span(index, index + rr), k) = c.slice(k)(arma::span(0, rr), rr);
            index = index + rr + 1;
        }
        index = 0;
    }
    return output;
}

inline
arma::mat cinrowmat (const arma::cube& c)
{
    const unsigned int r = c.n_rows;
    const unsigned int K = c.n_slices;
    arma::mat output(K, r * (r + 1) / 2);
    unsigned int index = 0;
    for (unsigned int k = 0; k < K; ++k) {
        for (unsigned int rr = 0; rr < r; ++rr) {
            output(k, arma::span(index, index + rr)) = arma::trans(c.slice(k)(arma::span(0, rr), rr));
            index = index + rr + 1;
        }
        index = 0;
    }
    return output;
}

inline
arma::mat qinmatr (const arma::rowvec& v) 
{
    const unsigned int s = v.n_elem;
    const unsigned int r = -0.5 + std::sqrt( 0.25 + 2 * s );
    arma::mat tmp(r, r);
    unsigned int index = 0;
    for ( unsigned int rr = 0; rr < r; ++rr) {
        tmp(arma::span(0, rr), rr) = arma::trans(v.cols(index, index + rr));
        tmp(rr, arma::span(0, rr)) = v.cols(index, index + rr);
        index = index + rr + 1;
    }
    return tmp;
}
#endif /* __FINMIX_MINCOL_H__ */



