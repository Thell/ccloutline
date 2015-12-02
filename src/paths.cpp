// paths.cpp: Cartesian coordinates of eEularian paths enclosing
//            connected components.
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

#include "paths.hpp"

namespace ccloutline {

cc_nested_sz_vec paths(cc_nested_sz_vec cc_cycles, const size_t& dim) {

  size_t a(0), b(0), offa(0), offb(0);
  for (auto& cycles : cc_cycles) {
    for (auto& cycle : cycles) {
      offa = cycle.size();
      offb = 2 * offa;
      cycle.resize(cycle.size() * 3);
      for( size_t i = 0; i < offa; i++ ) {
        b = cycle[i] / dim;
        a = cycle[i] % dim;
        std::swap(cycle[i + offa], a);
        std::swap(cycle[i + offb], b);
      }
    }
  }

  return cc_cycles;
}

} // namespace ccloutline

//' @title Connected component border cartesian paths.
//'
//' @description Connected component grid graph cartesian coordinate
//' paths for each eulerian cycle (outline) for each component.
//'
//' @param m A two dimensional numeric matrix.
//' @param reverse Begin at the minimal column major indexed vertex.
//'
//' @return Nested list of lists of arrays for all outlines of each
//' component such that \code{ccl_paths(mat)[[label]][[id]]}
//' yields an array in the form of
//' \code{[ v_1, y_1, x_1, v_2, y_2, x_2, ... , v_1, y_1, x_1]}.
//'
//' Each [[label]] index consists of at least 1 path and if the
//' component has holes there is an addition  path for each hole.
//'
//' @details There is a performance hit with `reverse`.
//'
//' @name ccl_paths
//' @seealso ccl_cycles
//' @export
// [[Rcpp::export]]
Rcpp::List ccl_paths(const arma::mat& m, const bool reverse = false) {
  const auto b(ccloutline::borders(m));
  const auto l(ccloutline::labels(b));
  const auto c(ccloutline::cycles(b, l, reverse));
  return Rcpp::wrap(ccloutline::paths(c, m.n_rows + 1));
}
