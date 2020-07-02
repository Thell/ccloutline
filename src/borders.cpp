// borders.cpp: 8-way connected border states .
//
// Copyright (C) 2015 Thell Fowler
//
// This file is part of ccloutline.
//
// ccloutline is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ccloutline is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ccloutline.  If not, see <http://www.gnu.org/licenses/>.

#include "borders.hpp"

namespace {

void validate_input(const arma::mat& m) {
  if (m.n_cols < 2 || m.n_rows < 2) {
    Rcpp::stop("Input matrix have >=2 rows and columns!");
  }
}

} // anonymous namespace

namespace ccloutline {

border_mat borders(const arma::mat& m) {
  using namespace arma;

  validate_input(m);

  // Column major border connection point values
  const border_t N(1), NW(2), W(4), SW(8), S(16), SE(32), E(64), NE(128);

  // Subview matrix dims.
  const uword r = m.n_rows - 1;
  const uword c = m.n_cols - 1;

  // vertical, horizontal, diagonal and antidiagonal equality.
  typedef border_mat bm;
  const bm v = conv_to<bm>::from(m.head_rows(r) == m.tail_rows(r));
  const bm h = conv_to<bm>::from(m.head_cols(c) == m.tail_cols(c));
  const bm d = conv_to<bm>::from(m(0, 0, size(r, c)) == m(1, 1, size(r, c)));
  const bm a = conv_to<bm>::from(m(1, 0, size(r, c)) == m(0, 1, size(r, c)));

  bm b = bm(size(m));
  b.fill(255);
  b.tail_rows(r) -= N * v;
  b.head_rows(r) -= S * v;
  b.tail_cols(c) -= W * h;
  b.head_cols(c) -= E * h;
  b(1, 1, size(r, c)) -= NW * d;
  b(0, 0, size(r, c)) -= SE * d;
  b(1, 0, size(r, c)) -= NE * a;
  b(0, 1, size(r, c)) -= SW * a;

  return b;
}

} // namespace ccloutline

//' @title 8-way connectedness border states.
//'
//' @description 8-way neighborhood connected component border
//' (edge/vertice) states.
//'
//' @param m A two dimensional numeric matrix.
//'
//' @return Numeric matrix of the same dimensions as mat with
//' values indicating 8-way connected border states.
//'
//' @details A cell is considered to have 8 distinct borders indexed
//' in counter clockwise order with values of 2^c(0:7) beginning
//' with the North edge and ending with the Northeast vertex. A cell
//' fully enclosed by similar neighbors has no borders ( a 0 value )
//' whereas when fully enclosed by dissimilar neighbors is fully
//' bordered ( a 255 value ).
//'
//' @name ccl_borders
//' @export
// [[Rcpp::export]]
 Rcpp::NumericMatrix ccl_borders(const arma::mat& m) {
  auto b(ccloutline::borders(m));
  arma::mat b_conv = arma::conv_to<arma::mat>::from(b);
  return Rcpp::wrap(b_conv);
}
