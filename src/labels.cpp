// labels.cpp: Connected components labelling and membership counts.
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

#include "labels.hpp"

#include <array>
#include <functional>

namespace ccloutline {

// Given a matrix of border states representing 8-way connected closed
// borders as set bits, non-connected = 255 and fully connected = 0,
// check the first four bits for openness to connected component.
// Use column major offsets from current index for ccl[ N, NW, W, SW ].
Labels labels(const border_mat& border_states) {
  using namespace arma;

  auto dim(border_states.n_rows);
  std::array<size_t, 4> ccl = {1, dim + 1, dim, dim - 1};

  sz_vec parents;
  parents.reserve(border_states.n_elem);
  std::vector<std::reference_wrapper<size_t>> labels;
  labels.reserve(border_states.n_elem);

  size_t index(0), label(0);
  border_t state(0);
  auto ccl_ref = std::ref(label);

  for (const auto& borders : border_states) {
    state = (~borders) & 15;
    if (state > 0) {
      ccl_ref = labels.at(index - ccl.at(__builtin_ffsl(state) - 1));
    } else {
      parents.emplace_back(parents.size() + 1);
      ccl_ref = parents.back();
    }
    labels.emplace_back(ccl_ref);

    if (state > 8 && state < 12) {
      rem_merge(parents, labels[index], labels.at(index - ccl[3]));
    }
    ++index;
  }

  // Flatten parent labels to a sequence.
  index = 0;
  label = 0;
  for (auto& p : parents) {
    p = p < ++index ? parents[p - 1] : ++label;
  }

  // The label matrix
  index = 0;
  label_mat ccl_umat = label_mat(size(border_states));
  ccl_umat.imbue([&]() { return labels[index++]; });

  // Component membership counts (count the refs to each parent).
  sz_vec member_counts(label);
  sz_vec p_indexes(parents.begin(), parents.end());
  for (auto& p : parents)
    p = 0;
  for (auto& p : labels)
    p = ++p;
  for (index = 0; index < parents.size(); ++index)
    member_counts[p_indexes[index] - 1] += parents[index];

  return std::make_pair(ccl_umat, member_counts);
}

} // namespace ccloutline

//' @title 8-way connected component labelling.
//'
//' @description 8-way neighborhood connected component labelling of
//' all given values.
//'
//' @param m A two dimensional numeric matrix.
//'
//' @return Numeric matrix of the same dimensions as mat with
//'   values indicating component identification label where values
//'   increase in column major order.
//'
//' @name ccl_labels
//' @export
// [[Rcpp::export]]
Rcpp::List ccl_labels(const arma::mat& m) {
  using namespace Rcpp;
  const auto b(ccloutline::borders(m));
  auto l(ccloutline::labels(b));
  arma::mat l_conv = arma::conv_to<arma::mat>::from(l.first);
  NumericVector l_sizes = wrap(l.second);
  return List::create(_["labels"]=wrap(l_conv),  _["sizes"]=l_sizes);
}
