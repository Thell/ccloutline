/*
 * See /test/profile/README.md
 *
 */

#include "ccloutline.hpp"

#ifdef CCL_GPROF
#include <gperftools/profiler.h>
#endif  // CCL_GPROF

// [[Rcpp::export]]
Rcpp::List ccl_profile( arma::mat input, bool reverse=false ) {

#ifdef CCL_GPROF
  ProfilerEnable();
  auto out = ccl_outlines(input, reverse);
  ProfilerDisable();
  return out;
#else
  return Rcpp::List::create(Rcpp::Named("NotEnabled")=0);
#endif

}

