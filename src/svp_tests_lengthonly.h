#ifndef SVPCHANGE_SVP_TESTS_LENGTHONLY_H
#define SVPCHANGE_SVP_TESTS_LENGTHONLY_H

#include "gaussian_focus_state.h"
#include "svp_test_base.h"

#include <cmath>
#include <cstddef>

// Length-only multiscale threshold (the SVP_validity_rules.tex proposal):
// the single-anchored statistic f_Gauss(y_{s:t}) = max_cost() compared
// against a threshold depending only on the candidate's own current length
// m = t - s (not on the split point), following the Dumbgen-Spokoiny/SMUCE
// multiscale correction sqrt(2 log(e*n/m)):
//
//   max_cost()/sigma2 <= [q + sqrt(2 * log(e * n / m))]^2
//
// where n is the TOTAL series length (fixed at construction) and q is the
// free calibration constant (passed as "gamma"). Dispatched with "both"
// pruning as a heuristic shortcut to this raw statistic, exactly like every
// other both-family method in this package: no F_both aggregate is
// computed. Unlike the original SVP_validity_rules.tex draft, every length
// m is scanned (no declared length subset), and a single shared q is used
// (no per-length gamma_m optimization).
class GaussianMeanLengthOnly : public TestBase
{
public:
  GaussianMeanLengthOnly(double sigma2, double n)
    : sigma2_(sigma2), n_(n), length_(0) {}

  void update(double y) override
  {
    focus_.update(y);
    ++length_;
  }

  std::size_t length() const { return length_; }

  double threshold(double q) const
  {
    const double m = static_cast<double>(length_);
    const double inner = q + std::sqrt(2.0 * std::log(kE * n_ / m));
    return inner * inner;
  }

  double statistic() const override
  {
    return focus_.max_cost() / sigma2_;
  }

  bool passes(double gamma) const override
  {
    return statistic() <= threshold(gamma);
  }

private:
  static constexpr double kE = 2.718281828459045235360287471352662498;

  double sigma2_;
  double n_;
  std::size_t length_;
  changepoint::GaussianFocusState focus_;
};

#endif
