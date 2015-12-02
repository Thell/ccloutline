// ccloutline.cpp: Cartesian coordinates of eulerian paths enclosing
//                 connected components.
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

#include "ccloutline.hpp"

//' @title Connected component border cartesian paths.
//'
//' @description Connected component grid graph cartesian coordinate
//' paths for each eulerian cycle (outline) for each component.
//'
//' @param mat A two dimensional numeric matrix.
//' @param min_first Begin at the minimal column major indexed vertex.
//'
//' @return Nested list of lists of matrices for all outlines of each
//' input label's connected components.
//'
//' Each [[label]] index consists of at least 1 path and if the
//' component has holes there is an additional path for each hole.
//'
//' @details There is a performance hit with `min_first`.
//'
//' @name ccl_outlines
//' @seealso ccl_borders, ccl_labels, ccl_vertices, ccl_cycles, ccl_paths
//' @export
// [[Rcpp::export]]
Rcpp::List ccl_outlines(const arma::mat& m, const bool min_first = true) {
  using namespace Rcpp;

  const auto b(ccloutline::borders(m));
  const auto l(ccloutline::labels(b));
  const auto c(ccloutline::cycles(b, l, min_first));
  const auto p(ccloutline::paths(c, m.n_rows + 1));

  arma::uvec li = arma::find_unique(l.first);

  List ccl_paths(p.size());
  auto j(0);
  for (const auto& ccl : p) {
    List paths(ccl.size());
    auto i(0);
    for (const auto& path : ccl) {
      arma::mat tmp = arma::conv_to<arma::mat>::from(path);
      tmp.reshape(path.size() / 3, 3);
      tmp.col(0).fill(i + 1);
      tmp.swap_cols(1, 2);
      tmp.insert_cols(0, 1);
      tmp.col(0).fill(m.at(li.at(j)));
      paths[i] = wrap(tmp);
      colnames(paths[i]) = CharacterVector::create("label", "ccl", "x", "y");
      i++;
    }
    ccl_paths[j] = paths;
    j++;
  }

  return ccl_paths;
}
