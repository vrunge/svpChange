#include "svp_tests_gaussian.h"
#include <Rcpp.h>
#include <vector>

// Internal C++ backend shared by valid_FOCUS() and valid_FOCUS_last().
// It remains outside an anonymous namespace because Rcpp exposes it through
// the generated .Call interface.
//
// [[Rcpp::export(name = ".focus_valid_cpp")]]
bool focus_valid_cpp(std::vector<double> data,
                     double gamma,
                     bool check_all_prefixes = true)
{
  GaussianMean test;
  for (double value : data) {
    test.update(value);
    if (check_all_prefixes && test.statistic() >= gamma) return false;
  }
  return test.statistic() < gamma;
}
