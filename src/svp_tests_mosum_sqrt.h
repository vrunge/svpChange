#ifndef SVPCHANGE_SVP_TESTS_MOSUM_SQRT_H
#define SVPCHANGE_SVP_TESTS_MOSUM_SQRT_H

#include "svp_test_base.h"

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

// MOSUM (moving-sum) validity test with a multiscale, length-adaptive
// penalty, replacing the FOCUS/CUSUM maximum used by every other test in
// this package with a bandwidth-G moving-sum statistic:
//
//   T_k(G) = (1/sqrt(2G)) * ( sum_{i=k+1}^{k+G} y_i - sum_{i=k-G+1}^{k} y_i )
//            / sigma,                        G <= k <= m - G,
//
// m = candidate's own current length. Unlike FOCUS, this scans local
// windows of FIXED width G rather than maximizing over every split jointly
// with the two segment means, so it is not itself a likelihood-ratio
// statistic for a single change -- it is the classical MOSUM scan of
// Husková and Slabý (2001, Kybernetika 37) for the i.i.d. multiple-change
// hypothesis.
//
// Multiscale penalty: the bandwidth G = G(m) = max(2, floor(sqrt(m)))
// grows with the candidate's own length m (so every candidate length is
// tested at an appropriately-scaled resolution, m/G -> infinity as
// m -> infinity, satisfying the regularity condition below), and the
// threshold on max_k T_k(G)^2 is the exact asymptotic boundary of
// Eichinger and Kirch (2018, Bernoulli 24(1), Theorem 2.1): for bandwidths
// G = G(n) with n/G -> infinity (and a mild rate condition relating G to
// n^{2/(2+nu)} log n, automatically satisfied here for Gaussian data),
//
//   a(n/G) * max_{G<=k<=n-G} |T_k(G)| - b(n/G)  -->  Gumbel(0,1) in
//   distribution,  P(Gumbel <= x) = exp(-2 exp(-x)),
//   a(x) = sqrt(2 log x),
//   b(x) = 2 log x + (1/2) log(log x) + log(3/2) - (1/2) log(pi).
//
// Substituting the candidate's own length m for n gives a length-adaptive
// (multiscale) threshold: solving the Gumbel CDF for its (1-alpha)
// quantile q_alpha = -log(-(1/2) log(1-alpha)) and inverting a(.)/b(.)
// gives threshold(m) = [(b(m/G) + q_alpha) / a(m/G)]^2 on the squared
// scale. statistic() returns max_k T_k(G)^2 - threshold(m) (the same
// "shift by an internal formula" convention as both_calibrated), so
// passes(gamma) with gamma = 0 reproduces the exact asymptotic boundary,
// and gamma is the single additive constant calibrated empirically (as
// for every other method in this document) to correct the inevitable
// finite-sample gap from the n -> infinity theory. Dispatched with "both"
// pruning as the same heuristic shortcut on this raw statistic used by
// every other both-family method here.
class GaussianMeanMosumSqrt : public TestBase
{
public:
  explicit GaussianMeanMosumSqrt(double sigma2 = 1.0, double alpha = 0.01)
    : sigma2_(sigma2), alpha_(alpha), prefix_(1, 0.0) {}

  void update(double y) override
  {
    prefix_.push_back(prefix_.back() + y);
  }

  std::size_t length() const { return prefix_.size() - 1; }

  static std::size_t bandwidth(std::size_t m)
  {
    const double g = std::floor(std::sqrt(static_cast<double>(m)));
    return std::max<std::size_t>(static_cast<std::size_t>(g), 2);
  }

  double raw_statistic(std::size_t m, std::size_t G) const
  {
    const double Gd = static_cast<double>(G);
    double best = -std::numeric_limits<double>::infinity();
    for (std::size_t k = G; k <= m - G; ++k) {
      const double right_sum = prefix_[k + G] - prefix_[k];
      const double left_sum  = prefix_[k] - prefix_[k - G];
      const double Tk = (right_sum - left_sum) / std::sqrt(2.0 * Gd * sigma2_);
      best = std::max(best, Tk * Tk);
    }
    return best;
  }

  double threshold(std::size_t m, std::size_t G) const
  {
    const double x = static_cast<double>(m) / static_cast<double>(G);
    const double logx = std::log(x);
    const double a = std::sqrt(2.0 * logx);
    const double b = 2.0 * logx + 0.5 * std::log(logx) +
      std::log(1.5) - 0.5 * std::log(M_PI);
    const double q = -std::log(-0.5 * std::log(1.0 - alpha_));
    const double thr = (b + q) / a;
    return thr * thr;
  }

  double statistic() const override
  {
    const std::size_t m = length();
    const std::size_t G = bandwidth(m);
    if (m < 2 * G) return -std::numeric_limits<double>::infinity();
    return raw_statistic(m, G) - threshold(m, G);
  }

  bool passes(double gamma) const override
  {
    return statistic() <= gamma;
  }

private:
  double sigma2_;
  double alpha_;
  std::vector<double> prefix_;
};

#endif
