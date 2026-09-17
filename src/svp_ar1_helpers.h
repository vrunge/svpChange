#ifndef SVPCHANGE_SVP_AR1_HELPERS_H
#define SVPCHANGE_SVP_AR1_HELPERS_H

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

// Shared implementation details for the AR(1) interfaces.
namespace svp_detail {

// Compute a median in-place. The input vector may be reordered by
// std::nth_element; this is used only while estimating AR(1) correlation.
inline double median_in_place(std::vector<double>& values)
{
  const std::size_t middle = values.size() / 2;
  std::nth_element(values.begin(), values.begin() + middle, values.end());
  const double upper = values[middle];
  if (values.size() % 2 == 1) return upper;
  const double lower = *std::max_element(values.begin(), values.begin() + middle);
  return 0.5 * (lower + upper);
}

// Estimate the AR(1) coefficient from the ratio of robust lag-two and lag-one
// difference scales.
inline double robust_ar1_rho(const std::vector<double>& data)
{
  if (data.size() < 3) {
    Rcpp::stop("AR(1) estimation requires at least 3 observations");
  }

  std::vector<double> lag1(data.size() - 1);
  std::vector<double> lag2(data.size() - 2);
  for (std::size_t i = 0; i < lag1.size(); ++i) {
    lag1[i] = std::fabs(data[i + 1] - data[i]);
  }
  for (std::size_t i = 0; i < lag2.size(); ++i) {
    lag2[i] = std::fabs(data[i + 2] - data[i]);
  }

  const double med1 = median_in_place(lag1);
  const double med2 = median_in_place(lag2);
  if (!std::isfinite(med1) || !std::isfinite(med2) || med1 <= 0.0) {
    Rcpp::stop(
      "Cannot estimate AR(1) correlation from constant or non-finite data"
    );
  }

  const double estimate = (med2 * med2) / (med1 * med1) - 1.0;
  // The model assumes stationarity. Clamping also prevents numerical
  // instability for short or highly contaminated series.
  return std::max(-0.999, std::min(0.999, estimate));
}

} // namespace svp_detail

#endif
