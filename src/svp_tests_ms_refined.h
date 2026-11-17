#ifndef SVPCHANGE_SVP_TESTS_MS_REFINED_H
#define SVPCHANGE_SVP_TESTS_MS_REFINED_H

#include "gaussian_focus_state.h"
#include "svp_test_base.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

// Single-change statistic on one interval. Under the SVP both rule, every
// contained interval must pass. The full length n is fixed for all intervals.
class GaussianMeanMSFamily : public TestBase {
public:
  GaussianMeanMSFamily(double sigma2, std::size_t n, double lambda,
                       double kappa)
    : sigma_(std::sqrt(sigma2)), n_(n), lambda_(lambda), kappa_(kappa),
      center_(0.0), seen_(false), prefix_(1, 0.0) {}

  void update(double y) override {
    if (!seen_) { center_ = y; seen_ = true; }
    const double centered = y - center_;
    focus_.update(centered);
    prefix_.push_back(prefix_.back() + centered);
  }

  std::size_t length() const { return prefix_.size() - 1; }

  double penalty(std::size_t m, std::size_t r) const {
    return std::sqrt(2.0 * (1.0 + lambda_ * std::log(
      static_cast<double>(m) / r) + kappa_ * std::log(
      static_cast<double>(n_) / m)));
  }

  double upper_bound() const {
    const std::size_t m = length();
    if (m < 2) return -std::numeric_limits<double>::infinity();
    return std::sqrt(focus_.max_cost()) / sigma_ - penalty(m, m / 2);
  }

  double statistic() const override {
    const std::size_t m = length();
    if (m < 2) return -std::numeric_limits<double>::infinity();
    double best = -std::numeric_limits<double>::infinity();
    const double total = prefix_[m];
    for (std::size_t left = 1; left < m; ++left) {
      const std::size_t right = m - left;
      const double z = std::abs(prefix_[left] / left -
        (total - prefix_[left]) / right) *
        std::sqrt(static_cast<double>(left) * right / m) / sigma_;
      best = std::max(best, z - penalty(m, std::min(left, right)));
    }
    return best;
  }

  bool passes(double q) const override {
    if (length() < 2) return true;
    const double margin = 64.0 * std::numeric_limits<double>::epsilon() *
      (1.0 + std::abs(q));
    return upper_bound() < q - margin || statistic() <= q;
  }

private:
  double sigma_;
  std::size_t n_;
  double lambda_;
  double kappa_;
  double center_;
  bool seen_;
  std::vector<double> prefix_;
  changepoint::GaussianFocusState focus_;
};

#endif
