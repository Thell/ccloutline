// edges.cpp: Connected component edge edges.
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

#include "edges.hpp"

#include <array>

namespace ccloutline {

// Edge vertices for each connected component label.
// [label][from,to,from,to,from,to,...]
nested_sz_vec edges(const border_mat& border_states, const Labels& ccl_stats) {

  auto ccl_count(ccl_stats.second.size());
  nested_sz_vec ccl_edges(ccl_count);
  for (size_t i = 0; i < ccl_count; ++i) {
    ccl_edges[i].reserve(8*ccl_stats.second[i]);
  }

  const auto NEdge(1), WEdge(4), SEdge(16), EEdge(64), Edge(85);
  size_t index(0), row(0), column(0);

  std::array<size_t, 4> vertex;
  auto update_vertex = [&]() {
    if (row == border_states.n_rows) {
      row = 0;
      ++column;
    }
    vertex[1] = index + column;
    vertex[2] = vertex[1] + 1;
    vertex[0] = vertex[2] + border_states.n_rows;
    vertex[3] = vertex[0] + 1;
  };

  for (const auto& state : border_states) {

    if (state & Edge) {

      update_vertex();
      auto& edges = ccl_edges[ccl_stats.first(index) - 1];

      if (state & NEdge) {
        edges.push_back(vertex[0]);
        edges.push_back(vertex[1]);
      }
      if (state & WEdge) {
        edges.push_back(vertex[1]);
        edges.push_back(vertex[2]);
      }
      if (state & SEdge) {
        edges.push_back(vertex[2]);
        edges.push_back(vertex[3]);
      }
      if (state & EEdge) {
        edges.push_back(vertex[3]);
        edges.push_back(vertex[0]);
      }
    }

    ++index;
    ++row;
  }

  return ccl_edges;
}

} // namespace ccloutline

//' @title Connected component border grid graph edges.
//'
//' @description Connected component directed edge vertex
//' connections.
//'
//' @param m A two dimensional numeric matrix.
//'
//' @return List of arrays for all directed edges for each connected
//' component such that \code{ccl_edges(mat)[[label]]} yields an
//' array in the form of
//' $[ from_1, to_1, from_2, to_2, ... , from_n, to_n]$.
//'
//' The first two sets of the returned array correspond to the
//' starting index (in column major order) of the component
//' within \code{m}.
//'
//' @name ccl_edges
//' @export
// [[Rcpp::export]]
Rcpp::List ccl_edges(const arma::mat& m) {
  const auto b(ccloutline::borders(m));
  const auto l(ccloutline::labels(b));
  return Rcpp::wrap(ccloutline::edges(b, l));
}
