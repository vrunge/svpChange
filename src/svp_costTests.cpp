#include "tests.h"
#include <Rcpp.h>
#include <functional>
#include <vector>
#include <limits>
#include <algorithm>
#include <iterator> // for std::make_move_iterator

using namespace Rcpp;

// [[Rcpp::export]]
List SVP_costTEsts(std::vector<double> data,
                  double gamma,
                  std::string test,
                  double quantile = 0.01)
{
  size_t n = data.size();

  NumericMatrix R(n + 1, 3); // Q, K, s
  R(0, 0) = 0.0;
  R(0, 1) = 0.0;
  R(0, 2) = 0.0;

  std::vector<size_t> nb(n); // nb of candidates examined at each t

  std::vector<double> S1(n + 1, 0);
  std::vector<double> S2(n + 1, 0);
  for (size_t i = 0; i < n; i++)
  {
    S1[i + 1] = S1[i] + data[i];
    S2[i + 1] = S2[i] + data[i] * data[i];
  }

  // factory for tests
  std::function<std::unique_ptr<TestBase>()> newTest;
  if (test == "gaussian_mean")
  {
    newTest = []() { return std::make_unique<GaussianMean>();};
  }
  else if (test == "gamma_rate")
  {
    newTest = []() { return std::make_unique<GammaRate>(); };
  }
  else if (test == "gaussian_variance")
  {
    newTest = []() { return std::make_unique<GaussianVariance>(); };
  }
  else if (test == "quantile")
  {
    newTest = [quantile]() { return std::make_unique<QuantileCostExact>(quantile); };
  }
  else if (test == "Chi2Cost")
  {
    newTest = []() { return std::make_unique<Chi2Cost>(); };
  }
  else
  {
    stop("Unknown test type");
  }





  std::vector<size_t> INDEX = {0};  // Active candidates (starting points s)
  std::vector<size_t> valid_INDEX;  // Candidates to keep for next iteration

  std::vector<std::unique_ptr<TestBase>> TESTS;  // Test instances for each candidate
  TESTS.push_back(newTest());

  // track how far each test has been updated (last index applied)
  std::vector<size_t> LAST_UPDATES;  // Last t value each test was updated to
  LAST_UPDATES.push_back(0);

  double best_Q;
  size_t best_K;
  size_t best_s = 0;
  size_t s;
  bool valid;
  double candidate_Q;
  size_t candidate_K;

    // PRUNING BRANCH: Keep only valid candidates and use lexicographic order
    for (size_t t = 1; t < n + 1; ++t)
    {
      // Rcpp::Rcout << "\n=== ITERATION t=" << t << " === Active candidates: " << INDEX.size() << std::endl;

      nb[t - 1] = INDEX.size();  // Record number of active candidates at this t

      best_Q = std::numeric_limits<double>::infinity();  // Best cost found
      best_K = std::numeric_limits<size_t>::max();  // Best K (lexicographic: primary criterion)

      valid_INDEX.clear();
      valid_INDEX.reserve(INDEX.size());  // Pre-allocate to avoid repeated allocations
      std::vector<std::unique_ptr<TestBase>> valid_TESTS;
      valid_TESTS.reserve(INDEX.size());
      std::vector<size_t> valid_LAST_UPDATES;
      valid_LAST_UPDATES.reserve(INDEX.size());

      for (size_t k = 0; k < INDEX.size(); ++k)
      {
        s = INDEX[k];  // Current candidate starting point
        auto& test_instance = TESTS[k];
        size_t last_up = LAST_UPDATES[k];

        // PRE-CHECK: Compute candidate_K using lexicographic order
        // K = number of segments = R(s,1) + 1. If K > best_K, skip remaining work
        candidate_K = static_cast<size_t>(R(s, 1)) + 1;

        // LEXICOGRAPHIC PRUNING: If candidate_K exceeds best_K found so far,
        // we can skip updating and evaluating this candidate AND all remaining ones
        // (because INDEX is kept in increasing order of R(s,1), so candidate_K will only increase)
        if (candidate_K > best_K) {
          // Rcpp::Rcout << "  [k=" << k << "] Lex-skip: candidate_K=" << candidate_K << " > best_K=" << best_K << " | Batch-appending k=" << k << "..." << (INDEX.size()-1) << std::endl;
          // IMPORTANT: do NOT permanently discard the current candidate or the remaining ones.
          // Append the current and all remaining candidates to the kept lists (without updating them now),
          // then break out of the loop. They will be evaluated (caught-up) in future t.
          valid_INDEX.insert(valid_INDEX.end(), INDEX.begin() + k, INDEX.end());

          valid_TESTS.insert(valid_TESTS.end(),
                             std::make_move_iterator(TESTS.begin() + k),
                             std::make_move_iterator(TESTS.end()));

          valid_LAST_UPDATES.insert(valid_LAST_UPDATES.end(),
                                    LAST_UPDATES.begin() + k,
                                    LAST_UPDATES.end());
          break;  // Exit loop, remaining candidates will be re-evaluated in future iterations
        }

        // CATCH-UP: Update test statistic from last_up to current t
        // (so test_instance->statistic() reflects the [s+1, t] segment)
        if (last_up < t) {
          // Rcpp::Rcout << "  [k=" << k << "] Catch-up from t=" << last_up << " to t=" << t << std::endl;
          for (size_t u = last_up + 1; u <= t; ++u) {
            test_instance->update(data[u - 1]);
          }
          last_up = t;
        }

        // EVALUATE CANDIDATE [s, t]
        candidate_Q = R(s, 0) + (S2[t] - S2[s]) - (S1[t] - S1[s]) * (S1[t] - S1[s]) / (t - s);
        candidate_K = static_cast<size_t>(R(s, 1)) + 1;  // K = # of segments

        // VALIDITY TEST: Check if segment [s+1, t] passes the statistical test
        valid = test_instance->statistic() < gamma;
        // Rcpp::Rcout << "  [k=" << k << "] s=" << s << ": valid=" << valid << " stat=" << test_instance->statistic() << " K=" << candidate_K << " Q=" << candidate_Q << std::endl;

        if (valid)
        {
          // KEEP: Valid candidate is retained for future iterations
          valid_INDEX.push_back(s);
          valid_TESTS.push_back(std::move(test_instance));
          valid_LAST_UPDATES.push_back(last_up);

          // LEXICOGRAPHIC UPDATE: Update best only if (K, Q) is lexicographically smaller
          // Primary: minimize K (# segments)
          // Secondary: minimize Q (cost) if K is tied
          if (candidate_K < best_K || (candidate_K == best_K && candidate_Q < best_Q))
          {
            // Rcpp::Rcout << "    -> NEW BEST: K=" << candidate_K << " (prev=" << best_K << "), Q=" << candidate_Q << " (prev=" << best_Q << ")" << std::endl;
            best_Q = candidate_Q;
            best_K = candidate_K;
            best_s = s;  // Best starting point for segment ending at t
          }
        }
        // else: PRUNED (invalid candidate discarded, does not contribute to best_*)
      }


      // STORE OPTIMAL SOLUTION AT t
      // Rcpp::Rcout << "  OPTIMAL[t=" << t << "]: K=" << best_K << " Q=" << best_Q << " (best_s=" << best_s << ")" << std::endl;
      R(t, 0) = best_Q;  // Optimal cost
      R(t, 1) = best_K;  // Optimal number of segments
      R(t, 2) = best_s;  // Optimal starting point of last segment

      // PRUNE: Keep only valid candidates for next iteration
      // This maintains the pruning strategy
      INDEX.swap(valid_INDEX);
      TESTS = std::move(valid_TESTS);
      LAST_UPDATES = std::move(valid_LAST_UPDATES);

      // ADD NEW CANDIDATE: s = t (potential starting point for future segments)
      TESTS.push_back(newTest());
      LAST_UPDATES.push_back(t);  // New test is up-to-date to index t (no obs applied yet)
      INDEX.push_back(t);
      // Rcpp::Rcout << "  Added candidate s=" << t << " | Total candidates: " << INDEX.size() << std::endl;
    }

  // BACKTRACKING: Reconstruct changepoints using optimal starting points
  // R(t, 2) = best_s = optimal starting point of segment ending at t
  std::vector<size_t> changepoints;
  changepoints.reserve(n / 2);  // Reserve space (typical case: O(log n) changepoints, but allocate safely)
  size_t i = n;
  // Rcpp::Rcout << "\n=== BACKTRACKING ==="  << std::endl;
  while (i > 0)
  {
    // Rcpp::Rcout << "  t=" << i << " -> best_s=" << R(i, 2) << std::endl;
    changepoints.push_back(i);
    i = R(i, 2);  // Move to start of previous segment
  }

  std::reverse(changepoints.begin(), changepoints.end());
  std::reverse(INDEX.begin(), INDEX.end());

  return List::create(
    _["changepoints"] = changepoints,
    _["lastIndexSet"] = INDEX,
    _["nb"] = nb,
    _["costQ"] = NULL,
    _["R"] = R(Range(1, R.nrow() - 1), Range(0, R.ncol() - 1))
  );
}
