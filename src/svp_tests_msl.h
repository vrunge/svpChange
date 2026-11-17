#ifndef SVPCHANGE_SVP_TESTS_MSL_H
#define SVPCHANGE_SVP_TESTS_MSL_H

#include "gaussian_focus_state.h"
#include "svp_test_base.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

// Exact single-change MS-l statistic for the current interval. The SVP
// "both" pruning rule combines these interval tests into hereditary validity:
// a segment survives only when all of its contained intervals pass.
class GaussianMeanMSL : public TestBase
{
public:
  explicit GaussianMeanMSL(double sigma2 = 1.0)
    : sigma_(std::sqrt(sigma2)), center_(0.0), seen_(false), prefix_(1, 0.0) {}

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

  // Every split pays at least the balanced-split penalty. The FOCUS maximum
  // is therefore a safe upper bound on the penalized maximum and lets most
  // intervals pass without scanning all splits.
  double upper_bound() const
  {
    const std::size_t m = length();
    if (m < 2) return -std::numeric_limits<double>::infinity();
    const double minimum_penalty = penalty(m, m / 2);
    return std::sqrt(focus_.max_cost()) / sigma_ - minimum_penalty;
  }

  double statistic() const override
  {
    const std::size_t m = length();
    if (m < 2) return -std::numeric_limits<double>::infinity();

    double best = -std::numeric_limits<double>::infinity();
    const double total = prefix_[m];
    const double size = static_cast<double>(m);
    for (std::size_t left = 1; left < m; ++left) {
      const std::size_t right = m - left;
      const double z = std::abs(prefix_[left] / left -
        (total - prefix_[left]) / right) *
        std::sqrt(static_cast<double>(left) * right / size) / sigma_;
      best = std::max(best, z - penalty(m, std::min(left, right)));
    }
    return best;
  }

  bool passes(double q) const override
  {
    if (length() < 2) return true;
    const double margin = 64.0 * std::numeric_limits<double>::epsilon() *
      (1.0 + std::abs(q));
    if (upper_bound() < q - margin) return true;
    return statistic() <= q;
  }

private:
  static double penalty(std::size_t length, std::size_t shorter_side)
  {
    return std::sqrt(2.0 * (1.0 + std::log(
      static_cast<double>(length) / shorter_side)));
  }

  double sigma_;
  double center_;
  bool seen_;
  std::vector<double> prefix_;
  changepoint::GaussianFocusState focus_;
};

#endif
