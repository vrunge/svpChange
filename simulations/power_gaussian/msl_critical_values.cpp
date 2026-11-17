#include <Rcpp.h>
#include "svp_tests_msl.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

// Full-series maximum over every interval and every internal split. Under the
// hereditary both rule, this is the smallest q allowing one valid segment.
// [[Rcpp::export]]
double gaussian_msl_critical_value_cpp(const Rcpp::NumericVector& data,
                                       double sigma2 = 1.0)
{
  const std::size_t n = data.size();
  if (n < 2 || !std::isfinite(sigma2) || sigma2 <= 0.0)
    Rcpp::stop("data must have at least two points and sigma2 must be positive");
  if (!std::all_of(data.begin(), data.end(),
                   [](double value) { return std::isfinite(value); }))
    Rcpp::stop("data must contain only finite values");

  std::vector<GaussianMeanMSL> intervals;
  intervals.reserve(n);
  double best = -std::numeric_limits<double>::infinity();

  for (std::size_t t = 0; t < n; ++t) {
    intervals.emplace_back(sigma2);
    for (std::size_t start = 0; start <= t; ++start) {
      GaussianMeanMSL& interval = intervals[start];
      interval.update(data[t]);
      if (interval.length() < 2) continue;
      const double margin = 64.0 * std::numeric_limits<double>::epsilon() *
        (1.0 + std::abs(best));
      if (interval.upper_bound() + margin >= best)
        best = std::max(best, interval.statistic());
    }
  }
  return best;
}
