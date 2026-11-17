#ifndef SVPCHANGE_SVP_VERZELEN_H
#define SVPCHANGE_SVP_VERZELEN_H

#include "gaussian_focus_state.h"
#include "svp_internal.h"
#include <deque>

namespace svp_detail {

// All indices are boundaries: (a,b,c) compares data[a:b] to data[b:c]
// in C++ indexing. N is fixed for the entire series, including subsegments.
class VerzelenScanner {
public:
  VerzelenScanner(const std::vector<double>& data, double sigma2, int N)
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
      Rcpp::stop("The data range is too large for the Verzelen statistic.");
    // Conservative roundoff allowance for prefix sums and FOCUS arithmetic.
    // If it overflows, the bound is disabled and every split is scanned.
    const double n = static_cast<double>(data.size());
    guard_ = 128.0 * std::numeric_limits<double>::epsilon() * n * n *
      (1.0 + energy) / sigma2_;
  }

  double extend(std::size_t a, std::size_t c, double inherited,
                const changepoint::GaussianFocusState& focus,
                bool accelerated = true) const
  {
    const std::size_t m = c - a;
    if (m < 2) return inherited;
    // The smallest penalty occurs at the most balanced split. FOCUS supplies
    // an upper bound on C^2; its active splits alone do NOT maximize our score.
    const double min_penalty = 2.0 *
      (log_N_ + logs_[m] - logs_[m / 2] - logs_[m - m / 2]);
    if (accelerated) {
      const double upper = focus.max_cost() / (4.0 * sigma2_) - min_penalty;
      if (std::isfinite(upper) && upper + guard_ < inherited) return inherited;
    }
    double best = inherited;
    for (std::size_t b = a + 1; b < c; ++b) {
      const double l = static_cast<double>(b - a);
      const double r = static_cast<double>(c - b);
      const double diff = (sums_[b] - sums_[a]) / l -
        (sums_[c] - sums_[b]) / r;
      const double z = diff / std::sqrt(sigma2_);
      const double score = (l * r / static_cast<double>(m)) * z * z / 4.0 -
        2.0 * (log_N_ + logs_[m] - logs_[b - a] - logs_[c - b]);
      if (!std::isfinite(score))
        Rcpp::stop("The data range is too large for the Verzelen statistic.");
      best = std::max(best, score);
    }
    return best;
  }

private:
  std::vector<double> sums_, logs_;
  double log_N_, sigma2_, guard_;
};

struct VerzelenCandidate {
  std::size_t boundary;
  double maximum;
  changepoint::GaussianFocusState focus;
  explicit VerzelenCandidate(std::size_t s)
    : boundary(s), maximum(-std::numeric_limits<double>::infinity()) {}
};

// Exact full-series maximum for calibration and internal reference checks.
inline double verzelen_statistic(const std::vector<double>& data,
                                 double sigma2, int N, bool accelerated)
{
  VerzelenScanner scanner(data, sigma2, N);
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

// Standalone SVP runner. Containment makes validity hereditary: once a
// candidate fails, it and every older boundary can be removed permanently.
inline Rcpp::List svp_verzelen_impl(const std::vector<double>& data,
                                   double gamma, double sigma2,
                                   const SegmentCost& segment_cost)
{
  const std::size_t n = data.size();
  VerzelenScanner scanner(data, sigma2, static_cast<int>(n));
  std::deque<VerzelenCandidate> candidates;
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
      // V_s(t) = max(V_s(t-1), V_{s+1}(t), H(s,t)).
      candidate.maximum = scanner.extend(s, t,
        std::max(candidate.maximum, suffix_max), candidate.focus);
      suffix_max = candidate.maximum;
      if (candidate.maximum > gamma) { remove = j; break; }
    }
    while (remove > 0) { candidates.pop_front(); --remove; }
    std::size_t best_K = n + 1, best_s = 0;
    double best_Q = std::numeric_limits<double>::infinity();
    // Ascending boundaries preserve the existing tie-breaking convention.
    for (const auto& candidate : candidates) {
      const auto s = candidate.boundary;
      const double q = Q[s] + segment_cost(s, t);
      if (is_better_candidate(K[s] + 1, q, best_K, best_Q)) {
        best_K = K[s] + 1; best_Q = q; best_s = s;
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
