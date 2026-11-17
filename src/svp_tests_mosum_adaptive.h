#ifndef SVPCHANGE_SVP_TESTS_MOSUM_ADAPTIVE_H
#define SVPCHANGE_SVP_TESTS_MOSUM_ADAPTIVE_H

#include "svp_test_base.h"

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

// Length-adaptive version of both_mosum. The raw statistic is the same
// single, fixed-midpoint two-sample test,
//
//   Delta(s, tau_mid, t), tau_mid = s + floor(m/2), m = candidate's own
//   current length,
//
// which, unlike f_Gauss, is EXACTLY chi^2_1-distributed at every fixed m
// (no internal maximum over a split is taken). What needs correcting is
// not the per-length marginal distribution (which never changes), but the
// REPEATED-TESTING multiplicity from re-evaluating the same candidate at
// every length it survives under "both" pruning.
//
// The correlation of the raw signed statistic Z_mid(m) across two lengths
// m < m' depends only on the ratio r = m'/m, and is COMPACTLY SUPPORTED:
// Corr(Z_mid(m), Z_mid(m')) = (2-r)/sqrt(r) for 1<=r<2, and exactly 0 for
// r>=2 (verified analytically from the Brownian-motion limit, and by
// simulation). After the log-time change u = log(m), the process
// {Z_mid(e^u)} is therefore STATIONARY with finite-range dependence --
// unlike f_Gauss's own running maximum, which is non-stationary and
// requires a single-log sqrt(2*log(m)) Darling-Erdos-type boundary. For a
// (locally non-degenerate) STATIONARY process, the classical theory for
// controlling the running maximum over repeated tests as the sample size
// grows is the "law of the iterated logarithm" repeated-significance-test
// framework (Darling and Erdos 1956, general boundary-crossing theory;
// Darling and Robbins 1967, 1968, repeated/sequential significance
// tests), which gives a much GENTLER, DOUBLE-LOG threshold growth than
// the single-log corrections used by every f_Gauss-based method in this
// document.
//
// This motivates the parametric family
//
//   gamma(m) = [q + sqrt(2*log(log(m/m0) + e))]^2, m0 = kMinLength = 20,
//
// (the "+e" inside the inner log regularizes it to exactly 1 at m=m0,
// avoiding log(0)). This is the same kind of "theoretically-motivated
// shape, empirically calibrated constant" construction as
// both_lengthonly (Dumbgen-Spokoiny shape) and both_mosum_sqrt (Eichinger-Kirch
// shape) elsewhere in this document: the exact extremal constant for this
// specific compactly-supported correlation structure is not derived here
// in closed form (unlike both_mosum_sqrt's Eichinger-Kirch boundary, which is
// exact); q is found empirically by simulation instead.
class GaussianMeanMosumAdaptive : public TestBase
{
public:
  explicit GaussianMeanMosumAdaptive(double sigma2 = 1.0)
    : sigma2_(sigma2), prefix_(1, 0.0) {}

  void update(double y) override
  {
    prefix_.push_back(prefix_.back() + y);
  }

  std::size_t length() const { return prefix_.size() - 1; }

  double raw_statistic() const
  {
    const std::size_t m = length();
    const std::size_t mid = m / 2;
    const double total = prefix_[m];
    const double left_sum = prefix_[mid];
    const double right_sum = total - left_sum;
    const double left_n = static_cast<double>(mid);
    const double right_n = static_cast<double>(m - mid);
    const double diff = left_sum / left_n - right_sum / right_n;
    return diff * diff * (left_n * right_n / static_cast<double>(m)) / sigma2_;
  }

  double threshold(double q) const
  {
    const std::size_t m = length();
    const double inner_log = std::log(
      static_cast<double>(m) / static_cast<double>(kMinLength));
    const double loglog = std::log(inner_log + M_E);
    const double t = q + std::sqrt(2.0 * loglog);
    return t * t;
  }

  double statistic() const override
  {
    if (length() < kMinLength) return -std::numeric_limits<double>::infinity();
    return raw_statistic();
  }

  bool passes(double gamma) const override
  {
    if (length() < kMinLength) return true;
    return raw_statistic() <= threshold(gamma);
  }

  static constexpr std::size_t kMinLength = 20;

private:
  double sigma2_;
  std::vector<double> prefix_;
};

#endif
