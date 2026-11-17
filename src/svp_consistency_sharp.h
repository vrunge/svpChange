#ifndef SVPCHANGE_SVP_CONSISTENCY_SHARP_H
#define SVPCHANGE_SVP_CONSISTENCY_SHARP_H

#include "gaussian_focus_state.h"
#include "svp_internal.h"
#include <deque>

namespace svp_detail {

// Implements the multiscale "both" validity test of
// simulations/consistency/SVP_consistency_sharp.tex, Section 6
// ("Multiscale validity with both subsegment tests"):
//
//   T_n(Y;s,t) = max_{s<=a<u<b<=t} { |Z_u(a,b)| - p_n(u-a,b-u) },
//   p_n(l,r)   = sqrt(2 * log(n*(l+r)/(l*r))),
//
// with n the TOTAL series length, fixed for the whole run (never the
// current candidate's own length -- the paper's nesting proof requires
// this, see its "Hereditary validity, implementation" subsection).
// A segment is valid iff T_n <= q. Structurally this is the same
// every-sub-interval scan as the existing "verzelen" test (same C(a,b,c)
// contrast, same containment/nesting argument for exact "both" pruning,
// Proposition "Stability for both pruning rules" in the paper), but on
// the |Z| (linear) scale with penalty p_n(l,r) instead of the Verzelen
// scan's squared, log(N*m/(l*r))-penalized statistic.
class ConsistencySharpScanner {
public:
  ConsistencySharpScanner(const std::vector<double>& data, double sigma2, int N)
    : sums_(data.size() + 1, 0.0), logs_(data.size() + 1, 0.0),
      log_N_(std::log(static_cast<double>(N))), sigma2_(sigma2), guard_(0.0)
  {
    double energy = 0.0;
    for (std::size_t i = 0; i < data.size(); ++i) {
      const double x = data[i] - data[0];
      sums_[i + 1] = sums_[i] + x;
      energy += x * x;
      logs_[i + 1] = std::log(static_cast<double>(i + 1));
    }
    if (!std::isfinite(energy) || !std::isfinite(sums_.back()))
      Rcpp::stop("The data range is too large for the consistency-sharp statistic.");
    const double n = static_cast<double>(data.size());
    // Conservative roundoff allowance, propagated through the square root
    // from the same order-of-magnitude bound used by the Verzelen scanner.
    guard_ = std::sqrt(128.0 * std::numeric_limits<double>::epsilon() * n * n *
      (1.0 + energy) / sigma2_);
  }

  double extend(std::size_t a, std::size_t c, double inherited,
                const changepoint::GaussianFocusState& focus,
                bool accelerated = true) const
  {
    const std::size_t m = c - a;
    if (m < 2) return inherited;
    // p_n(l,r) is minimized, for fixed m = l + r, at the most balanced
    // split (l*r maximal). |Z_b(a,c)| <= sqrt(max_cost()/sigma2) for every
    // b (an exact identity, since max_cost() = max_b Z_b(a,c)^2). Their
    // difference therefore bounds every term of the exact scan below.
    const double min_penalty = std::sqrt(2.0 *
      (log_N_ + logs_[m] - logs_[m / 2] - logs_[m - m / 2]));
    if (accelerated) {
      const double upper = std::sqrt(focus.max_cost() / sigma2_) - min_penalty;
      if (std::isfinite(upper) && upper + guard_ < inherited) return inherited;
    }
    double best = inherited;
    for (std::size_t b = a + 1; b < c; ++b) {
      const double l = static_cast<double>(b - a);
      const double r = static_cast<double>(c - b);
      const double diff = (sums_[b] - sums_[a]) / l -
        (sums_[c] - sums_[b]) / r;
      const double z = diff / std::sqrt(sigma2_);
      const double abs_z = std::sqrt(l * r / static_cast<double>(m)) * std::fabs(z);
      const double penalty = std::sqrt(2.0 *
        (log_N_ + logs_[m] - logs_[b - a] - logs_[c - b]));
      const double score = abs_z - penalty;
      if (!std::isfinite(score))
        Rcpp::stop("The data range is too large for the consistency-sharp statistic.");
      best = std::max(best, score);
    }
    return best;
  }

private:
  std::vector<double> sums_, logs_;
  double log_N_, sigma2_, guard_;
};

struct ConsistencySharpCandidate {
  std::size_t boundary;
  double maximum;
  changepoint::GaussianFocusState focus;
  explicit ConsistencySharpCandidate(std::size_t s)
    : boundary(s), maximum(-std::numeric_limits<double>::infinity()) {}
};

// Exact full-series maximum, for calibration (simulating M_n = T_n(sigma*eps;0,n)
// per the paper's calibration recipe: "simulate the maximum M_n itself").
inline double consistency_sharp_statistic(const std::vector<double>& data,
                                          double sigma2, int N, bool accelerated)
{
  ConsistencySharpScanner scanner(data, sigma2, N);
  double best = -std::numeric_limits<double>::infinity();
  for (std::size_t a = 0; a < data.size(); ++a) {
    Rcpp::checkUserInterrupt();
    changepoint::GaussianFocusState focus;
    for (std::size_t c = a + 1; c <= data.size(); ++c) {
      focus.update(data[c - 1] - data[a]);
      best = scanner.extend(a, c, best, focus, accelerated);
    }
  }
  return best;
}

// Standalone SVP runner. Containment makes validity hereditary here exactly
// as for "verzelen" (Proposition "Stability for both pruning rules"): once
// a candidate fails, it and every older boundary can be removed permanently.
inline Rcpp::List svp_consistency_sharp_impl(const std::vector<double>& data,
                                             double q, double sigma2,
                                             const SegmentCost& segment_cost)
{
  const std::size_t n = data.size();
  ConsistencySharpScanner scanner(data, sigma2, static_cast<int>(n));
  std::deque<ConsistencySharpCandidate> candidates;
  candidates.emplace_back(0);
  std::vector<double> Q(n + 1, 0.0);
  std::vector<std::size_t> K(n + 1, 0), previous(n + 1, 0), nb(n, 0);
  for (std::size_t t = 1; t <= n; ++t) {
    Rcpp::checkUserInterrupt();
    nb[t - 1] = candidates.size();
    double suffix_max = -std::numeric_limits<double>::infinity();
    std::size_t remove = 0;
    for (std::size_t j = candidates.size(); j > 0; --j) {
      auto& candidate = candidates[j - 1];
      const std::size_t s = candidate.boundary;
      candidate.focus.update(data[t - 1] - data[s]);
      candidate.maximum = scanner.extend(s, t,
        std::max(candidate.maximum, suffix_max), candidate.focus);
      suffix_max = candidate.maximum;
      if (candidate.maximum > q) { remove = j; break; }
    }
    while (remove > 0) { candidates.pop_front(); --remove; }
    std::size_t best_K = n + 1, best_s = 0;
    double best_Q = std::numeric_limits<double>::infinity();
    for (const auto& candidate : candidates) {
      const auto s = candidate.boundary;
      const double val = Q[s] + segment_cost(s, t);
      if (is_better_candidate(K[s] + 1, val, best_K, best_Q)) {
        best_K = K[s] + 1; best_Q = val; best_s = s;
      }
    }
    Q[t] = best_Q; K[t] = best_K; previous[t] = best_s;
    candidates.emplace_back(t);
  }
  std::vector<std::size_t> changepoints, last_index_set;
  for (std::size_t t = n; t > 0; t = previous[t]) changepoints.push_back(t);
  std::reverse(changepoints.begin(), changepoints.end());
  for (auto it = candidates.rbegin(); it != candidates.rend(); ++it)
    last_index_set.push_back(it->boundary);
  return Rcpp::List::create(
    Rcpp::Named("changepoints") = changepoints,
    Rcpp::Named("lastIndexSet") = last_index_set,
    Rcpp::Named("nb") = nb,
    Rcpp::Named("costQ") = R_NilValue,
    Rcpp::Named("R") = build_R_matrix(Q, K, previous));
}

} // namespace svp_detail
#endif
