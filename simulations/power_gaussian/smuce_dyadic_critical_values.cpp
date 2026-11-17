#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>

// Smallest nonnegative dyadic-length SMUCE threshold making the full series
// a valid one-segment fit under the "dyaLen" interval system (every starting
// position, but only dyadic lengths 1, 2, 4, ...). Mirrors
// smuce_critical_values.cpp exactly except the outer loop only visits
// dyadic lengths, giving O(n log n) instead of O(n^2) work. Used only to
// calibrate the Gaussian null.
// [[Rcpp::export]]
double smuce_dyadic_critical_value_cpp(const Rcpp::NumericVector& data,
                                       double sigma2 = 1.0) {
  const int n = data.size();
  if (n < 1 || !R_finite(sigma2) || sigma2 <= 0.0)
    Rcpp::stop("invalid arguments");

  const double offset = data[0];
  std::vector<double> prefix(n + 1, 0.0);
  for (int i = 0; i < n; ++i) {
    if (!R_finite(data[i])) Rcpp::stop("data must be finite");
    prefix[i + 1] = prefix[i] + data[i] - offset;
  }

  std::vector<double> intercept_down, intercept_up, slope;
  for (int length = 1; length <= n; length *= 2) {
    double max_sum = R_NegInf;
    double min_sum = R_PosInf;
    for (int start = 0; start + length <= n; ++start) {
      const double sum = prefix[start + length] - prefix[start];
      max_sum = std::max(max_sum, sum);
      min_sum = std::min(min_sum, sum);
    }
    const double root = std::sqrt(static_cast<double>(length) / sigma2);
    const double penalty = std::sqrt(
      2.0 * (1.0 + std::log(static_cast<double>(n) / length))
    );
    slope.push_back(root);
    intercept_down.push_back(max_sum / std::sqrt(length * sigma2) - penalty);
    intercept_up.push_back(-min_sum / std::sqrt(length * sigma2) - penalty);
  }
  const int m = static_cast<int>(slope.size());

  // The one-segment statistic is min_theta max(g(theta), h(theta)), where
  // g is the upper envelope of decreasing lines and h of increasing lines,
  // now built from dyadic lengths only. Their intersection is the
  // minimizer; each envelope evaluation is O(m) = O(log n).
  double low = *std::min_element(data.begin(), data.end()) - offset;
  double high = *std::max_element(data.begin(), data.end()) - offset;
  for (int iter = 0; iter < 64; ++iter) {
    const double theta = low + (high - low) / 2.0;
    double down = R_NegInf;
    double up = R_NegInf;
    for (int i = 0; i < m; ++i) {
      down = std::max(down, intercept_down[i] - slope[i] * theta);
      up = std::max(up, intercept_up[i] + slope[i] * theta);
    }
    if (down > up) low = theta;
    else high = theta;
  }

  const double theta = low + (high - low) / 2.0;
  double value = R_NegInf;
  for (int i = 0; i < m; ++i) {
    value = std::max(value, intercept_down[i] - slope[i] * theta);
    value = std::max(value, intercept_up[i] + slope[i] * theta);
  }
  return std::max(0.0, value);
}
