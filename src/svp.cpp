#include "tests.h"
#include <Rcpp.h>
#include <functional>
#include <vector>
#include <limits>
#include <algorithm>

using namespace Rcpp;

// [[Rcpp::export]]
List SVP(std::vector<double> data,
         double gamma,
         std::string test,
         bool prune_if_unvalid = true)
{
  size_t n = data.size();

  NumericMatrix R(n + 1, 3); // Q, K, s
  R(0, 0) = 0.0;
  R(0, 1) = 0.0;
  R(0, 2) = 0.0;

  std::vector<size_t> nb(n); // nb of candidates examined at each t

  // cumulative sums
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
  else
  {
    stop("Unknown test type");
  }

  std::vector<size_t> INDEX = {0};
  std::vector<size_t> valid_INDEX;

  std::vector<std::unique_ptr<TestBase>> TESTS;
  TESTS.push_back(newTest());

  // track how far each test has been updated (last index applied)
  std::vector<size_t> LAST_UPDATES;
  LAST_UPDATES.push_back(0);

  double best_Q;
  size_t best_K;
  size_t best_s = 0;
  size_t s;
  bool valid;
  double candidate_Q;
  size_t candidate_K;

  if (prune_if_unvalid == true)
  {
    for (size_t t = 1; t < n + 1; ++t)
    {
      nb[t - 1] = INDEX.size();

      best_Q = std::numeric_limits<double>::infinity();
      best_K = std::numeric_limits<size_t>::max();

      valid_INDEX.clear();
      std::vector<std::unique_ptr<TestBase>> valid_TESTS;
      std::vector<size_t> valid_LAST_UPDATES;

      for (size_t k = 0; k < INDEX.size(); ++k)
      {
        s = INDEX[k];
        auto& test_instance = TESTS[k];
        size_t last_up = LAST_UPDATES[k];

        // compute candidate_K and possibly skip (safe because instance is up-to-date)
        candidate_K = static_cast<size_t>(R(s, 1)) + 1;
        if (candidate_K > best_K) {
          continue; // no chance to improve lexicographic order
        }


        // catch up this instance up to current t (so statistic() is correct)
        if (last_up < t) {
          for (size_t u = last_up + 1; u <= t; ++u) {
            test_instance->update(data[u - 1]);
          }
          last_up = t;
        }

        candidate_Q = R(s, 0) + (S2[t] - S2[s]) - (S1[t] - S1[s]) * (S1[t] - S1[s]) / (t - s);
        candidate_K = static_cast<size_t>(R(s, 1)) + 1;

        // evaluate validity (using up-to-date stats)
        valid = test_instance->statistic() < gamma;

        if (valid)
        {
          valid_INDEX.push_back(s);
          valid_TESTS.push_back(std::move(test_instance));
          valid_LAST_UPDATES.push_back(last_up);



          if (candidate_K < best_K || (candidate_K == best_K && candidate_Q < best_Q))
          {
            best_Q = candidate_Q;
            best_K = candidate_K;
            best_s = s;
          }
        }
        // else: pruned (do not use invalid candidates to update best_*)
      }


      R(t, 0) = best_Q;
      R(t, 1) = best_K;
      R(t, 2) = best_s;

      // replace with kept valid lists and add new entry for s = t
      INDEX.swap(valid_INDEX);
      TESTS = std::move(valid_TESTS);
      LAST_UPDATES = std::move(valid_LAST_UPDATES);

      TESTS.push_back(newTest());
      LAST_UPDATES.push_back(t); // new test is up-to-date to index t (no obs applied)
      INDEX.push_back(t);
    }
  }
  else
  {
    // Non-pruning branch kept consistent with LAST_UPDATES maintenance
    for (size_t t = 1; t < n + 1; ++t)
    {
      nb[t - 1] = INDEX.size();

      best_Q = std::numeric_limits<double>::infinity();
      best_K = std::numeric_limits<size_t>::max();

      for (size_t k = 0; k < INDEX.size(); ++k)
      {
        s = INDEX[k];
        auto& test_instance = TESTS[k];

        size_t last_up = LAST_UPDATES[k];
        if (last_up < t) {
          for (size_t u = last_up + 1; u <= t; ++u) {
            test_instance->update(data[u - 1]);
          }
          LAST_UPDATES[k] = t;
        }

        valid = test_instance->statistic() < gamma;

        if (valid == true)
        {
          candidate_Q = R(s, 0) + (S2[t] - S2[s]) - (S1[t] - S1[s]) * (S1[t] - S1[s]) / (t - s);
          candidate_K = static_cast<size_t>(R(s, 1)) + 1;

          if (candidate_K < best_K || (candidate_K == best_K && candidate_Q < best_Q))
          {
            best_Q = candidate_Q;
            best_K = candidate_K;
            best_s = s;
          }
        }
      }
      R(t, 0) = best_Q;
      R(t, 1) = best_K;
      R(t, 2) = best_s;

      TESTS.push_back(newTest());
      LAST_UPDATES.push_back(t);
      INDEX.push_back(t);
    }
  }

  // BACKTRACKING
  std::vector<size_t> changepoints;
  size_t i = n;
  while (i > 0)
  {
    changepoints.push_back(i);
    i = R(i, 2);
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
