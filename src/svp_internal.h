#ifndef SVPCHANGE_SVP_INTERNAL_H
#define SVPCHANGE_SVP_INTERNAL_H

#include <Rcpp.h>
#include "svp_test_base.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <utility>
#include <vector>

// Private implementation of the SVP dynamic program. The templates are kept
// in this header so that svp.cpp can instantiate svp_impl() for each validity
// test used by the public dispatcher.
namespace svp_detail {

// Within-segment cost used to rank partitions that have the same number of
// segments. Boundaries follow the SVP convention: s is zero-based and t is
// one-based, so the pair (s, t) is the R segment data[(s + 1):t].
//
// Two cost series are supported. The Gaussian cost sums the squared deviations
// of data[s + 1:t] from the segment mean. The AR(1) cost does the same for the
// innovations z[u] = data[u] - rho * data[u - 1], which are the residuals of
// the AR(1) mean model: inside a segment with constant mean mu they have a
// constant mean (1 - rho) * mu, so their sum of squared deviations is the
// profiled conditional cost of that segment.
//
// The AR(1) cost skips z[s + 1], the innovation that straddles the boundary.
// That innovation has mean mu_new - rho * mu_old rather than
// (1 - rho) * mu_new, so keeping it would charge every segment for the jump
// that precedes it and pull boundaries away from the change. SVP compares
// costs only within a fixed number of segments K, and every partition into K
// segments skips exactly K - 1 innovations, so the comparison stays fair.
class SegmentCost
{
public:
  // Gaussian cost over data[s + 1:t].
  static SegmentCost gaussian(const std::vector<double>& data)
  {
    return SegmentCost(data, 0, 0.0);
  }

  // AR(1) innovation cost over z[s + 2], ..., z[t].
  static SegmentCost ar1(const std::vector<double>& data, double rho)
  {
    return SegmentCost(data, 1, rho);
  }

  double operator()(std::size_t s, std::size_t t) const
  {
    // A segment shorter than the offset contributes no cost term: it is a
    // single observation under the Gaussian cost, or a segment with no
    // interior innovation under the AR(1) cost.
    if (t < s + offset_ + 1) return 0.0;

    const std::size_t first = s + offset_;
    const double sum = S1_[t] - S1_[first];
    const double count = static_cast<double>(t - first);
    const double cost = (S2_[t] - S2_[first]) - sum * sum / count;

    if (cost <= 0.0) return 0.0;
    return cost;
  }

private:
  // offset is the number of leading observations of a segment that carry no
  // cost term; rho is used only when offset is one.
  SegmentCost(const std::vector<double>& data, std::size_t offset, double rho)
    : offset_(offset), S1_(data.size() + 1, 0.0), S2_(data.size() + 1, 0.0)
  {
    const std::size_t n = data.size();

    // Terms are indexed by their one-based endpoint u, and the first
    // "offset_" of them do not exist. term(u) is the Gaussian observation
    // data[u] or the AR(1) innovation data[u] - rho * data[u - 1].
    std::vector<double> term(n + 1, 0.0);
    for (std::size_t u = offset_ + 1; u <= n; ++u) {
      term[u] = offset_ == 0 ? data[u - 1] : data[u - 1] - rho * data[u - 2];
    }

    // Center only the cost calculation. This does not alter the validity test
    // or the data passed to it, but avoids cancellation in
    // sum(term^2) - sum(term)^2 / count.
    double center = 0.0;
    std::size_t seen = 0;
    for (std::size_t u = offset_ + 1; u <= n; ++u) {
      ++seen;
      center += (term[u] - center) / static_cast<double>(seen);
    }
    if (!std::isfinite(center)) {
      Rcpp::stop("Unable to compute a stable center for the segment costs.");
    }

    for (std::size_t u = offset_ + 1; u <= n; ++u) {
      const double centered = term[u] - center;
      if (!std::isfinite(centered)) {
        Rcpp::stop("The data range is too large for stable segment costs.");
      }
      S1_[u] = S1_[u - 1] + centered;
      S2_[u] = S2_[u - 1] + centered * centered;
      if (!std::isfinite(S1_[u]) || !std::isfinite(S2_[u])) {
        Rcpp::stop("The data range is too large for stable segment costs.");
      }
    }
  }

  std::size_t offset_;
  std::vector<double> S1_;
  std::vector<double> S2_;
};

// State associated with one possible previous boundary in the dynamic
// program. Each candidate owns an incremental validity-test object for its
// segment.
template <typename Test>
struct CandidateState
{
  std::size_t boundary;    // zero-based start boundary of the candidate segment
  Test test;               // incremental validity-test state
  std::size_t last_update; // last endpoint already processed by test
  bool valid;               // validity at the current endpoint

  template <typename... Args>
  CandidateState(std::size_t boundary_, Args&&... args)
    : boundary(boundary_),
      test(std::forward<Args>(args)...),
      last_update(boundary_),
      valid(true)
  {}

  CandidateState(const CandidateState&) = delete;
  CandidateState& operator=(const CandidateState&) = delete;
  CandidateState(CandidateState&&) = default;
  CandidateState& operator=(CandidateState&&) = default;
};

// Extend a candidate to endpoint t by feeding it only observations that have
// not been processed yet, then evaluate its validity statistic at t. Singleton
// segments are always valid, as in svp0().
template <typename Test>
void update_to_endpoint(CandidateState<Test>& candidate,
                        const std::vector<double>& data,
                        std::size_t t,
                        double gamma)
{
  // All built-in SVP tests are evaluated at the current endpoint. A previous
  // failure does not permanently invalidate a candidate, which is required
  // when subtests == "none".
  candidate.valid = true;

  for (std::size_t u = candidate.last_update + 1; u <= t; ++u) {
    candidate.test.update(data[u - 1]);

    // As in svp0(), a one-observation segment is always valid.
    const bool check_statistic =
      u - candidate.boundary > 1 &&
      u == t;
    if (check_statistic) {
      if (!candidate.test.passes(gamma)) {
        candidate.valid = false;
        candidate.last_update = t;
        return;
      }
    }
  }
  candidate.last_update = t;
}

// Reconstruct a deferred right-pruning candidate. Unlike update_to_endpoint(),
// this checks every missed endpoint because right pruning would have removed
// the candidate permanently at the first invalid endpoint. The replay is
// needed only when the valid one-segment shortcut stops being applicable.
template <typename Test, typename... Args>
void replay_right_candidate(CandidateState<Test>& candidate,
                            const std::vector<double>& data,
                            std::size_t t,
                            double gamma,
                            Args&... args)
{
  candidate.test = Test(args...);
  candidate.last_update = candidate.boundary;
  candidate.valid = true;

  for (std::size_t u = candidate.boundary + 1; u <= t; ++u) {
    candidate.test.update(data[u - 1]);
    candidate.last_update = u;

    if (u - candidate.boundary > 1 &&
        !candidate.test.passes(gamma)) {
      candidate.valid = false;
      return;
    }
  }
}

// Compare two valid partitions lexicographically: fewer segments are
// preferred first, and SSE is used only when the segment counts are equal.
inline bool is_better_candidate(std::size_t candidate_K,
                                double candidate_Q,
                                std::size_t best_K,
                                double best_Q)
{
  return candidate_K < best_K ||
    (candidate_K == best_K && candidate_Q < best_Q);
}

// Convert the internal dynamic-programming vectors into the three-column R
// matrix: cumulative cost Q, segment count K, and previous boundary.
inline Rcpp::NumericMatrix build_R_matrix(
    const std::vector<double>& Q,
    const std::vector<std::size_t>& K,
    const std::vector<std::size_t>& previous)
{
  const std::size_t n = Q.size() - 1;
  Rcpp::NumericMatrix R(n, 3);
  auto output = R.begin();

  // R matrices are stored column-major. Fill each complete column so that
  // writes are contiguous in memory; index zero is the DP initialization row.
  std::copy(Q.begin() + 1, Q.end(), output);
  std::transform(
    K.begin() + 1, K.end(), output + n,
    [](std::size_t value) {
      return static_cast<double>(value);
    }
  );
  std::transform(
    previous.begin() + 1, previous.end(), output + 2 * n,
    [](std::size_t value) {
      return static_cast<double>(value);
    }
  );

  return R;
}

// Run SVP for one concrete validity-test type. "cost" ranks partitions that
// have the same number of segments. Test-specific constructor arguments are
// forwarded through Args (for example quantile, rho, or sigma2). The public
// SVP() function validates arguments and selects the Test type and the cost
// before calling this implementation.
template <typename Test, typename... Args>
Rcpp::List svp_impl(const std::vector<double>& data,
                    double gamma,
                    const std::string& subtests,
                    const SegmentCost& cost,
                    Args&&... args)
{
  const std::size_t n = data.size();

  std::vector<double> Q(n + 1, std::numeric_limits<double>::infinity());
  std::vector<std::size_t> K(
    n + 1, std::numeric_limits<std::size_t>::max()
  );
  std::vector<std::size_t> previous(n + 1, 0);
  Q[0] = 0.0;
  K[0] = 0;

  std::vector<std::size_t> nb(n);

  using Candidate = CandidateState<Test>;
  std::vector<Candidate> candidates;
  std::vector<Candidate> retained;
  candidates.reserve(n + 1);
  retained.reserve(n + 1);
  candidates.emplace_back(0, std::forward<Args>(args)...);

  const bool no_pruning = subtests == "none";
  const bool prune_right_only = subtests == "right";

  // In the unpruned algorithm, candidate boundaries are never deleted. Group
  // them by their already-computed predecessor segment count K_s. At endpoint
  // t, the first group containing a valid candidate determines the smallest
  // possible K_t; only SSE ties within that group still need comparison.
  std::vector<std::vector<std::size_t>> candidates_by_K;
  if (no_pruning) {
    candidates_by_K.resize(n + 1);
    candidates_by_K[0].push_back(0);
  }

  // Dynamic-programming loop over endpoints t = 1, ..., n.
  for (std::size_t t = 1; t <= n; ++t) {
    nb[t - 1] = candidates.size();

    // With no candidate pruning, evaluate predecessor groups in increasing
    // K_s and stop after the first group containing a valid segment. All
    // candidates in that group are compared in increasing boundary order, so
    // the pre-existing strict SSE tie rule is unchanged.
    if (no_pruning) {
      double best_Q = std::numeric_limits<double>::infinity();
      std::size_t best_K = std::numeric_limits<std::size_t>::max();
      std::size_t best_s = 0;

      for (std::size_t predecessor_K = 0;
           predecessor_K < candidates_by_K.size(); ++predecessor_K) {
        bool found_valid = false;

        for (std::size_t boundary : candidates_by_K[predecessor_K]) {
          Candidate& candidate = candidates[boundary];
          update_to_endpoint(candidate, data, t, gamma);
          if (!candidate.valid) continue;

          const double candidate_Q = Q[boundary] + cost(boundary, t);
          if (!std::isfinite(candidate_Q)) continue;

          found_valid = true;
          if (candidate_Q < best_Q) {
            best_Q = candidate_Q;
            best_K = predecessor_K + 1;
            best_s = boundary;
          }
        }

        if (found_valid) break;
      }

      if (!std::isfinite(best_Q)) {
        Rcpp::stop("No valid finite-cost partition was found at the current endpoint.");
      }

      Q[t] = best_Q;
      K[t] = best_K;
      previous[t] = best_s;
      candidates.emplace_back(t, std::forward<Args>(args)...);
      candidates_by_K[best_K].push_back(t);
      continue;
    }

    // Equation (4): with right pruning, a valid candidate starting at zero
    // gives a one-segment partition, which is unbeatable in the primary
    // objective. Defer all later candidates until this candidate first fails.
    // When they are eventually processed, update_to_endpoint() catches their
    // test state up but evaluates validity only at the current endpoint.
    // This shortcut is deliberately restricted to "right": in "both", a
    // later invalid candidate can trigger left pruning and invalidate the
    // older one-segment candidate.
    if (prune_right_only && !candidates.empty() &&
        candidates.front().boundary == 0) {
      update_to_endpoint(candidates.front(), data, t, gamma);
      if (candidates.front().valid) {
        Q[t] = cost(0, t);
        K[t] = 1;
        previous[t] = 0;
        candidates.emplace_back(t, std::forward<Args>(args)...);
        continue;
      }

      // The shortcut has ended. Reconstruct the exact right-pruning history
      // for candidates that were deferred while the one-segment candidate was
      // valid. Candidates invalid at any missed endpoint are discarded.
      retained.clear();
      for (std::size_t i = 1; i < candidates.size(); ++i) {
        replay_right_candidate(candidates[i], data, t, gamma, args...);
        if (candidates[i].valid) {
          retained.push_back(std::move(candidates[i]));
        }
      }
      candidates.swap(retained);
    }

    double best_Q = std::numeric_limits<double>::infinity();
    std::size_t best_K = std::numeric_limits<std::size_t>::max();
    std::size_t best_s = 0;
    retained.clear();

    if (prune_right_only) {
      for (Candidate& candidate : candidates) {
        update_to_endpoint(candidate, data, t, gamma);
        if (!candidate.valid) continue;

        retained.push_back(std::move(candidate));
        const std::size_t boundary = retained.back().boundary;
        if (K[boundary] == std::numeric_limits<std::size_t>::max()) continue;

        const std::size_t candidate_K = K[boundary] + 1;

        // Valid survivors must be retained for future endpoints, but an SSE
        // cannot improve an already smaller segment count. Skip only that
        // arithmetic; validity and candidate-set diagnostics stay unchanged.
        if (best_K != std::numeric_limits<std::size_t>::max() &&
            candidate_K > best_K) {
          continue;
        }

        const double candidate_Q = Q[boundary] + cost(boundary, t);
        if (std::isfinite(candidate_Q) &&
            is_better_candidate(candidate_K, candidate_Q,
                                best_K, best_Q)) {
          best_Q = candidate_Q;
          best_K = candidate_K;
          best_s = boundary;
        }
      }
    } else {
      // "both" combines right pruning with the inclusive left rule. Because
      // boundaries are ordered increasingly, scan from newest to oldest and
      // stop at the first invalid candidate: it and every older candidate are
      // deleted, so their validity tests and costs cannot affect the result.
      std::size_t first_retained = 0;
      for (std::size_t i = candidates.size(); i-- > 0;) {
        update_to_endpoint(candidates[i], data, t, gamma);
        if (!candidates[i].valid) {
          first_retained = i + 1;
          break;
        }
      }

      // All surviving candidates are valid. Identify the smallest attainable
      // segment count before computing SSE, then compare SSE only inside that
      // class. The ascending pass preserves the strict-tie preference for the
      // smallest boundary used by the previous implementation.
      for (std::size_t i = first_retained; i < candidates.size(); ++i) {
        const std::size_t boundary = candidates[i].boundary;
        if (K[boundary] != std::numeric_limits<std::size_t>::max()) {
          best_K = std::min(best_K, K[boundary] + 1);
        }
      }
      for (std::size_t i = first_retained; i < candidates.size(); ++i) {
        const std::size_t boundary = candidates[i].boundary;
        if (K[boundary] == std::numeric_limits<std::size_t>::max() ||
            K[boundary] + 1 != best_K) {
          continue;
        }

        const double candidate_Q = Q[boundary] + cost(boundary, t);
        if (std::isfinite(candidate_Q) && candidate_Q < best_Q) {
          best_Q = candidate_Q;
          best_s = boundary;
        }
      }
      for (std::size_t i = first_retained; i < candidates.size(); ++i) {
        retained.push_back(std::move(candidates[i]));
      }
    }

    if (!std::isfinite(best_Q)) {
      Rcpp::stop("No valid finite-cost partition was found at the current endpoint.");
    }

    Q[t] = best_Q;
    K[t] = best_K;
    previous[t] = best_s;

    candidates.swap(retained);

    candidates.emplace_back(t, std::forward<Args>(args)...);
  }

  std::vector<std::size_t> changepoints;
  changepoints.reserve(n / 2 + 1);
  for (std::size_t i = n; i > 0;) {
    if (previous[i] >= i) {
      Rcpp::stop("Internal error while reconstructing the SVP partition.");
    }
    changepoints.push_back(i);
    i = previous[i];
  }
  std::reverse(changepoints.begin(), changepoints.end());

  std::vector<std::size_t> last_index_set;
  last_index_set.reserve(candidates.size());
  for (auto it = candidates.rbegin(); it != candidates.rend(); ++it) {
    last_index_set.push_back(it->boundary);
  }

  return Rcpp::List::create(
    Rcpp::Named("changepoints") = changepoints,
    Rcpp::Named("lastIndexSet") = last_index_set,
    Rcpp::Named("nb") = nb,
    Rcpp::Named("costQ") = R_NilValue,
    Rcpp::Named("R") = build_R_matrix(Q, K, previous)
  );
}

} // namespace svp_detail

#endif
