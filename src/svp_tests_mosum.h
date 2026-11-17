#ifndef SVPCHANGE_SVP_TESTS_MOSUM_H
#define SVPCHANGE_SVP_TESTS_MOSUM_H

#include "svp_test_base.h"

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

// MOSUM at the maximal bandwidth G = floor(m/2): the general MOSUM
// statistic T_k(G) (see GaussianMeanMosumSqrt in svp_tests_mosum_sqrt.h)
// compares the mean of G points to the right of k against G points to the
// left, maximized over every admissible k in [G, m-G]. Setting G to its
// largest possible value, G = floor(m/2), collapses that range to a
// SINGLE point, k = floor(m/2) -- the segment's own midpoint -- so there
// is nothing left to maximize over: this class is literally MOSUM's
// degenerate, maximal-bandwidth special case, not a separate construction.
// Instead of maximizing the Gaussian generalized-likelihood-ratio
// statistic f_Gauss(y_{s:t}) over every internal split tau (as every
// other test in this package does), it evaluates the split at that one
// fixed midpoint position tau = floor(m/2), m = t - s the candidate's
// current length:
//
//   Delta(s, tau_mid, t) = (ybar_right - ybar_left)^2 * (left*right/m) / sigma2,
//   tau_mid = s + floor(m/2).
//
// Segments shorter than kMinLength (20 observations) are always valid
// (never tested): below that length the fixed midpoint split leaves too
// few points on one side to be informative, so this is a hard floor, not
// an asymptotic one.
class GaussianMeanMosum : public TestBase
{
public:
  explicit GaussianMeanMosum(double sigma2 = 1.0)
    : sigma2_(sigma2), prefix_(1, 0.0) {}

  void update(double y) override
  {
    prefix_.push_back(prefix_.back() + y);
  }

  std::size_t length() const { return prefix_.size() - 1; }

  double statistic() const override
  {
    const std::size_t m = length();
    if (m < kMinLength) return -std::numeric_limits<double>::infinity();
    const std::size_t mid = m / 2;
    const double total = prefix_[m];
    const double left_sum = prefix_[mid];
    const double right_sum = total - left_sum;
    const double left_n = static_cast<double>(mid);
    const double right_n = static_cast<double>(m - mid);
    const double diff = left_sum / left_n - right_sum / right_n;
    return diff * diff * (left_n * right_n / static_cast<double>(m)) / sigma2_;
  }

  bool passes(double gamma) const override
  {
    return statistic() <= gamma;
  }

  static constexpr std::size_t kMinLength = 20;

private:
  double sigma2_;
  std::vector<double> prefix_;
};

#endif
