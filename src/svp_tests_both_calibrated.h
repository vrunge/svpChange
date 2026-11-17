#ifndef SVPCHANGE_SVP_TESTS_BOTH_CALIBRATED_H
#define SVPCHANGE_SVP_TESTS_BOTH_CALIBRATED_H

#include "gaussian_focus_state.h"
#include "svp_test_base.h"

#include <cmath>
#include <cstddef>
#include <limits>

// Deliberately mismatched calibration, kept to make the mismatch concrete.
//
// statistic() is f_FOCUS(y_{s:t}) = max_{s<tau<t} Z^2_{s,tau,t}: the SAME
// single-anchored statistic as plain GaussianMean (focus_.max_cost()) --
// start s and end t are the candidate's own fixed
// endpoints; only the internal split tau is maxed over. This is NOT
// max_{s<=u<v<=t} f_FOCUS(y_{u:v}) (also maxing over u and v): that double
// max over every sub-interval of [s,t] is a different, more expensive
// statistic this class does not compute.
//
// The threshold, however, is
//
//   gamma(l) = A * log(l) + B * log(log(l)) + C,   l = candidate's own length
//
// with (A, B, C) fit to the empirical 99% null quantile of
// max_{1<=u<v<=n} f_FOCUS(y_{u:v}) -- the double max over every sub-interval
// of the TOTAL series -- as a function of n. That
// quantile is the right calibration for the double-max statistic evaluated
// at length l, not for the single-anchored f_FOCUS(y_{s:t}) this class
// actually compares it to. The two statistics differ (the double max is
// never smaller), so gamma(l) is too large for what's being tested here:
// expected to be badly over-conservative under "none" pruning (empirically,
// close to 0% null rejection instead of 1%). Segments shorter than
// kMinLength are always valid.
//
// Because gamma(l) grows with l while max_cost() also grows with l, nothing
// guarantees their difference is monotone under widening, so "both"'s
// pruning shortcut (invalidate-on-failure) is not formally justified here.
// Dispatched with "both" pruning anyway at the caller's direction; results
// should be checked against the unpruned "none" rule if exactness is in
// question.
class GaussianMeanBothCalibrated : public TestBase
{
public:
  explicit GaussianMeanBothCalibrated(
      double sigma2 = 1.0,
      double A = kDefaultA, double B = kDefaultB, double C = kDefaultC)
    : sigma2_(sigma2), length_(0), A_(A), B_(B), C_(C) {}

  void update(double y) override
  {
    focus_.update(y);
    ++length_;
  }

  std::size_t length() const { return length_; }

  double penalty(std::size_t m) const
  {
    const double log_m = std::log(static_cast<double>(m));
    return A_ * log_m + B_ * std::log(log_m) + C_;
  }

  double statistic() const override
  {
    if (length_ < kMinLength) return -std::numeric_limits<double>::infinity();
    return focus_.max_cost() / sigma2_ - penalty(length_);
  }

  bool passes(double gamma) const override
  {
    return statistic() <= gamma;
  }

  static constexpr double kDefaultA = 1.682969;
  static constexpr double kDefaultB = 6.208626;
  static constexpr double kDefaultC = 4.939248;

private:
  static constexpr std::size_t kMinLength = 10;

  double sigma2_;
  std::size_t length_;
  double A_, B_, C_;
  changepoint::GaussianFocusState focus_;
};

#endif
