#include <Rcpp.h>
#include "gaussian_focus_state.h"
#include <algorithm>
#include <cstddef>
#include <deque>

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
    const double statistic = prefix_state.max_cost();
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
          both_critical, state.max_cost()
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

// Running trajectory of M_b = max_{u<tau<b} Z(u,tau,b)^2 for every prefix
// length b = 2, ..., n, computed in one pass (same O(n^2) total cost as
// gaussian_critical_values_cpp's single endpoint value, since it maintains
// one GaussianFocusState per start and reads its max_cost() at every step,
// but exposes every intermediate b instead of only b = n). Index i
// (0-based) of the returned vector holds M_{i+2}.
//
// M_b fixes the RIGHT endpoint at b and maxes only over the start u and
// internal split tau: this is F_left(y_{0..b}), NOT F_both(y_{0..b}).
// F_both(y_{0..b}) -- the statistic that actually needs controlling for
// "both" pruning, matching what gaussian_critical_values_cpp's own "$both"
// computes via its full double loop over every start AND every
// intermediate endpoint -- is F_both(y_{0..b}) = max_{b'<=b} F_left(y_{0..b'}),
// i.e. the running cummax of this function's own output over b. Callers
// that need F_both (e.g. calibrate_empirical_both.R) must apply cummax()
// themselves; this function intentionally returns the raw, un-cummax'd
// F_left ingredient.
//
// [[Rcpp::export]]
Rcpp::NumericVector gaussian_both_trajectory_cpp(
    const Rcpp::NumericVector& data)
{
  const std::size_t n = static_cast<std::size_t>(data.size());
  if (n < 2) return Rcpp::NumericVector(0);

  std::deque<changepoint::GaussianFocusState> candidates;
  Rcpp::NumericVector trajectory(n - 1);

  for (std::size_t t = 0; t < n; ++t) {
    candidates.emplace_back();
    for (auto& c : candidates) c.update(data[t]);
    double running_max = 0.0;
    for (std::size_t i = 0; i + 1 < candidates.size(); ++i) {
      running_max = std::max(running_max, candidates[i].max_cost());
    }
    if (t >= 1) trajectory[t - 1] = running_max; // b = t + 1 >= 2
  }
  return trajectory;
}
