#ifndef SVPCHANGE_SVP_TESTS_RMIN_H
#define SVPCHANGE_SVP_TESTS_RMIN_H

#include "gaussian_focus_state.h"
#include "svp_test_base.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

// Two-sample Gaussian CUSUM test restricted to splits whose shorter arm has
// at least r_min observations, and at least frac times the interval's
// current length. The SVP "both" pruning rule combines these interval tests
// into hereditary validity: a segment survives only when every contained
// interval passes. Unlike a multiscale penalty, a split that fails either
// floor is not down-weighted: it is dropped from the scan entirely and can
// never cause a candidate to fail. The absolute floor removes the
// short-scale tests that otherwise dominate the null maximum of a
// constant-threshold scan; the relative floor additionally removes
// lopsided splits of long intervals, whose shorter arm is short only
// relative to the interval. The cost is that a change cannot be resolved
// within r_eff(length) observations of either end of its segment, where
// r_eff(m) = max(r_min, ceil(frac * m)).
class GaussianMeanRMin : public TestBase
{
public:
  GaussianMeanRMin(double sigma2, std::size_t r_min, double frac = 0.0)
    : sigma2_(sigma2), r_min_(std::max<std::size_t>(r_min, 1)), frac_(frac),
      center_(0.0), seen_(false), prefix_(1, 0.0) {}

  void update(double y) override
  {
    if (!seen_) {
      center_ = y;
      seen_ = true;
    }
    const double centered = y - center_;
    focus_.update(centered);
    prefix_.push_back(prefix_.back() + centered);
  }

  std::size_t length() const { return prefix_.size() - 1; }

  std::size_t effective_floor(std::size_t m) const
  {
    if (frac_ <= 0.0) return r_min_;
    const std::size_t relative = static_cast<std::size_t>(
      std::ceil(frac_ * static_cast<double>(m)));
    return std::max(r_min_, relative);
  }

  bool has_eligible_split() const
  {
    const std::size_t m = length();
    return m >= 2 * effective_floor(m);
  }

  // The unconstrained FOCuS maximum bounds the r_min-restricted maximum
  // from above, because the restricted scan is a subset of all splits. Most
  // intervals pass without the exact O(m) scan below.
  double upper_bound() const
  {
    if (!has_eligible_split())
      return -std::numeric_limits<double>::infinity();
    return focus_.max_cost() / sigma2_;
  }

  double statistic() const override
  {
    const std::size_t m = length();
    const std::size_t floor = effective_floor(m);
    if (m < 2 * floor) return -std::numeric_limits<double>::infinity();

    double best = -std::numeric_limits<double>::infinity();
    const double total = prefix_[m];
    const double size = static_cast<double>(m);
    for (std::size_t left = floor; left <= m - floor; ++left) {
      const std::size_t right = m - left;
      const double diff = prefix_[left] / left -
        (total - prefix_[left]) / right;
      const double z2 = diff * diff *
        (static_cast<double>(left) * right / size) / sigma2_;
      best = std::max(best, z2);
    }
    return best;
  }

  bool passes(double gamma) const override
  {
    if (!has_eligible_split()) return true;
    const double margin = 64.0 * std::numeric_limits<double>::epsilon() *
      (1.0 + std::abs(gamma));
    if (upper_bound() < gamma - margin) return true;
    return statistic() <= gamma;
  }

private:
  double sigma2_;
  std::size_t r_min_;
  double frac_;
  double center_;
  bool seen_;
  std::vector<double> prefix_;
  changepoint::GaussianFocusState focus_;
};

#endif
