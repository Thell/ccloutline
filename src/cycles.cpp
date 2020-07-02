// cycles.cpp: Eulerian path vertices enclosing connected components.
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

#include "cycles.hpp"

namespace ccloutline {

size_t add_edges(const sz_vec& edges, Edge_Table& edge_table) {
  size_t edge_count(0), limit(edges.size());
  for (size_t i = 0; i < limit; i += 2, edge_count++) {
    auto& edge = edge_table[edges[i]];
    std::get<0>(edge) += 1;
    if (std::get<0>(edge) == 1) {
      std::get<1>(edge) = edges[i + 1];
    } else {
      std::get<2>(edge) = edges[i + 1];
    }
  }
  return edge_count;
};

size_t starting_vertex(size_t lbound, const Edge_Table& edge_table) {
  auto limit(edge_table.size());
  while (lbound < limit && std::get<0>(edge_table[lbound]) == 0) {
    lbound++;
  };
  return lbound;
};

size_t advance(incident_t& incident_to, vertex_t vertex, incident_vec& i_stack,
               sz_vec& v_stack, Edge_Table& edge_table) {
  v_stack.emplace_back(vertex);
  vertex = incident_to == 1 ? std::get<1>(edge_table[vertex])
                            : std::get<2>(edge_table[vertex]);
  --incident_to;
  i_stack.emplace_back(incident_to);
  return vertex;
};

size_t retreat(vertex_t vertex, incident_vec& i_stack, sz_vec& v_stack,
               sz_vec& cycle) {
  cycle.emplace_back(vertex);
  while (!v_stack.empty() && i_stack.back() == 0) {
    cycle.emplace_back(v_stack.back());
    v_stack.pop_back();
    i_stack.pop_back();
  }
  if (!v_stack.empty()) {
    vertex = v_stack.back();
    v_stack.pop_back();
    i_stack.pop_back();
  }
  return vertex;
};

size_t populate_cycle(size_t& edge_count, vertex_t vertex,
                      incident_vec& i_stack, sz_vec& v_stack, sz_vec& cycle,
                      Edge_Table& edge_table) {
  do {
    auto& incident_to = std::get<0>(edge_table[vertex]);
    if (incident_to == 0) {
      vertex = retreat(vertex, i_stack, v_stack, cycle);
    } else {
      vertex = advance(incident_to, vertex, i_stack, v_stack, edge_table);
      --edge_count;
    }
  } while (!v_stack.empty());
  return edge_count;
};

cc_nested_sz_vec cycles(const border_mat& borders, const Labels& labels,
                        const bool reverse = false) {

  const auto max_v((borders.n_cols + 1) * (borders.n_rows + 1));

  Edge_Table edge_table(max_v);

  sz_vec v_stack;
  v_stack.reserve(max_v);

  incident_vec i_stack;
  i_stack.reserve(max_v);

  cc_nested_sz_vec cc_cycles;
  cc_cycles.reserve(labels.second.size());

  nested_sz_vec cycles;

  sz_vec cycle;
  cycle.reserve(max_v);

  size_t cc_lbound(0), cycle_lbound(0), remaining_edges(0);

  for (const auto& cc_edges : edges(borders, labels)) {

    remaining_edges = add_edges(cc_edges, edge_table);
    cc_lbound = starting_vertex(cc_lbound, edge_table);
    cycle_lbound = cc_lbound;

    while (remaining_edges > 0) {
      cycle_lbound = starting_vertex(cycle_lbound, edge_table);
      remaining_edges = populate_cycle(remaining_edges, cycle_lbound, i_stack,
                                       v_stack, cycle, edge_table);
      if (reverse)
        std::reverse(cycle.begin(), cycle.end());
      cycles.emplace_back(cycle);
      cycle.clear();
    }

    cc_cycles.emplace_back(cycles);
    cycles.clear();
  }

  return cc_cycles;
}

} // namespace ccloutline

//' @title Connected component border vertice paths.
//'
//' @description Connected component grid graph vertice paths for
//' each eulerian cycle for each component.
//'
//' @param m A two dimensional numeric matrix.
//' @param reverse Begin at the minimal column major indexed vertex.
//'
//' @return Nested list of lists of arrays for all outlines of each
//' component such that \code{ccl_cycles(mat)[[label]][[id]]} yields
//' an array in the form of
//' \code{[ v_1, v_2, v_3, ... , v_1]}.
//'
//' Each [[label]] index consists of at least 1 path and if the
//' component has holes there is an additional path for each hole.
//'
//' @details There is a slight performance hit with `reverse` as this
//' forces a scan of the vertices for the minimal remaining valued
//' vertice and reorders the final results.
//'
//' @name ccl_cycles
//' @export
// [[Rcpp::export]]
Rcpp::List ccl_cycles(const arma::mat& m, bool reverse = false) {
  const auto b(ccloutline::borders(m));
  const auto l(ccloutline::labels(b));
  auto c(ccloutline::cycles(b, l, reverse));
  return Rcpp::wrap(c);
}
