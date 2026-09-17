#include <Rcpp.h>
#include "gaussian_focus_state.h"
#include <algorithm>
#include <cstddef>

// Return the null thresholds above which the complete series is a valid
// one-segment solution for each subtests mode. This simulation-only helper
// avoids fitting SVP repeatedly over a grid of constants.
//
// [[Rcpp::export]]
Rcpp::NumericVector gaussian_critical_values_cpp(
    const Rcpp::NumericVector& data)
{
  const std::size_t n = static_cast<std::size_t>(data.size());
  changepoint::GaussianFocusState prefix_state;
  double right_critical = 0.0;
  double none_critical = 0.0;

  for (std::size_t t = 0; t < n; ++t) {
    prefix_state.update(data[t]);
    const double statistic = 0.5 * prefix_state.max_cost();
    right_critical = std::max(right_critical, statistic);
    if (t + 1 == n) none_critical = statistic;
  }

  double both_critical = 0.0;
  for (std::size_t s = 0; s < n; ++s) {
    changepoint::GaussianFocusState state;
    for (std::size_t t = s; t < n; ++t) {
      state.update(data[t]);
      if (t > s) {
        both_critical = std::max(
          both_critical, 0.5 * state.max_cost()
        );
      }
    }
  }

  return Rcpp::NumericVector::create(
    Rcpp::Named("none") = none_critical,
    Rcpp::Named("right") = right_critical,
    Rcpp::Named("both") = both_critical
  );
}
