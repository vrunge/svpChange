#ifndef SVPCHANGE_SVP_TESTS_FOCUS_LIL_H
#define SVPCHANGE_SVP_TESTS_FOCUS_LIL_H

#include "gaussian_focus_state.h"
#include "svp_test_base.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

// FOCuS computes max_l Z_l^2 exactly. The length-only penalty is the maximum
// of the position penalty in Verzelen et al. (2023), Eq. (12), times L^2=4.
// Thus this statistic is bounded by their Proposition 2 statistic pointwise.
class GaussianMeanFocusLIL : public TestBase {
public:
  explicit GaussianMeanFocusLIL(double sigma2 = 1.0)
    : sigma2_(sigma2), center_(0.0), seen_(false), length_(0) {}

  void update(double y) override {
    if (!seen_) { center_ = y; seen_ = true; }
    focus_.update(y - center_);
    ++length_;
  }

  std::size_t length() const { return length_; }

  // max_{1<=l<m} max(min(l,m/l),min(m-l,m/(m-l))) is attained at a
  // split nearest sqrt(m), hence this exact two-point calculation.
  static double penalty(std::size_t m) {
    const std::size_t lower = static_cast<std::size_t>(
      std::sqrt(static_cast<double>(m)));
    const std::size_t upper = lower * lower == m ? lower : lower + 1;
    const double h = std::max(static_cast<double>(lower),
                              static_cast<double>(m) / upper);
    return 8.0 * std::log1p(std::log(h));
  }

  double statistic() const override {
    if (length_ < 2) return -std::numeric_limits<double>::infinity();
    return focus_.max_cost() / sigma2_ - penalty(length_);
  }

  bool passes(double q) const override {
    return length_ < 2 || statistic() <= q;
  }

private:
  double sigma2_;
  double center_;
  bool seen_;
  std::size_t length_;
  changepoint::GaussianFocusState focus_;
};

#endif
