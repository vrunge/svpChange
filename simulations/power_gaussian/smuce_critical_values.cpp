#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>

// Smallest nonnegative SMUCE threshold making the full series a valid
// one-segment fit. This is used only to calibrate the Gaussian null.
// [[Rcpp::export]]
double smuce_critical_value_cpp(const Rcpp::NumericVector& data,
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

  std::vector<double> intercept_down(n), intercept_up(n), slope(n);
  for (int length = 1; length <= n; ++length) {
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
    slope[length - 1] = root;
    intercept_down[length - 1] = max_sum / std::sqrt(length * sigma2) - penalty;
    intercept_up[length - 1] = -min_sum / std::sqrt(length * sigma2) - penalty;
  }

  // The one-segment statistic is min_theta max(g(theta), h(theta)), where
  // g is the upper envelope of decreasing lines and h of increasing lines.
  // Their intersection is the minimizer; each envelope evaluation is O(n).
  double low = *std::min_element(data.begin(), data.end()) - offset;
  double high = *std::max_element(data.begin(), data.end()) - offset;
  for (int iter = 0; iter < 64; ++iter) {
    const double theta = low + (high - low) / 2.0;
    double down = R_NegInf;
    double up = R_NegInf;
    for (int i = 0; i < n; ++i) {
      down = std::max(down, intercept_down[i] - slope[i] * theta);
      up = std::max(up, intercept_up[i] + slope[i] * theta);
    }
    if (down > up) low = theta;
    else high = theta;
  }

  const double theta = low + (high - low) / 2.0;
  double value = R_NegInf;
  for (int i = 0; i < n; ++i) {
    value = std::max(value, intercept_down[i] - slope[i] * theta);
    value = std::max(value, intercept_up[i] + slope[i] * theta);
  }
  return std::max(0.0, value);
}
