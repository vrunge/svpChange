#include "tests.h"
#include <Rcpp.h>
#include <vector>
#include <limits>
#include <algorithm>

using namespace Rcpp;

namespace {

inline double segment_cost(const std::vector<double>& S1,
                           const std::vector<double>& S2,
                           size_t s,
                           size_t t)
{
  const double sum = S1[t] - S1[s];
  return (S2[t] - S2[s]) - sum * sum / static_cast<double>(t - s);
}

NumericMatrix build_R_matrix(const std::vector<double>& Q,
                             const std::vector<size_t>& K,
                             const std::vector<size_t>& previous)
{
  const size_t n = Q.size() - 1;
  NumericMatrix R(n, 3);
  for (size_t t = 1; t <= n; ++t) {
    const size_t row = t - 1;
    R(row, 0) = Q[t];
    R(row, 1) = static_cast<double>(K[t]);
    R(row, 2) = static_cast<double>(previous[t]);
  }
  return R;
}

template <typename Test>
List svp_impl(const std::vector<double>& data,
              double gamma,
              bool prune_if_unvalid)
{
  const size_t n = data.size();

  std::vector<double> Q(n + 1, std::numeric_limits<double>::infinity());
  std::vector<size_t> K(n + 1, std::numeric_limits<size_t>::max());
  std::vector<size_t> previous(n + 1, 0);
  Q[0] = 0.0;
  K[0] = 0;

  std::vector<size_t> nb(n);

  std::vector<double> S1(n + 1, 0.0);
  std::vector<double> S2(n + 1, 0.0);
  for (size_t i = 0; i < n; ++i) {
    S1[i + 1] = S1[i] + data[i];
    S2[i + 1] = S2[i] + data[i] * data[i];
  }

  std::vector<size_t> index;
  std::vector<Test> tests;
  std::vector<size_t> last_updates;
  index.reserve(n + 1);
  tests.reserve(n + 1);
  last_updates.reserve(n + 1);
  index.push_back(0);
  tests.emplace_back();
  last_updates.push_back(0);

  for (size_t t = 1; t <= n; ++t) {
    const size_t m = index.size();
    nb[t - 1] = m;

    double best_Q = std::numeric_limits<double>::infinity();
    size_t best_K = std::numeric_limits<size_t>::max();
    size_t best_s = 0;

    if (prune_if_unvalid) {
      size_t write = 0;

      for (size_t k = 0; k < m; ++k) {
        const size_t s = index[k];
        const size_t candidate_K = K[s] + 1;

        if (candidate_K > best_K) {
          for (size_t j = k; j < m; ++j, ++write) {
            if (write != j) {
              index[write] = index[j];
              tests[write] = std::move(tests[j]);
              last_updates[write] = last_updates[j];
            }
          }
          break;
        }

        size_t last_up = last_updates[k];
        for (size_t u = last_up + 1; u <= t; ++u) {
          tests[k].update(data[u - 1]);
        }
        last_up = t;

        if (tests[k].statistic() < gamma) {
          const double candidate_Q = Q[s] + segment_cost(S1, S2, s, t);

          if (candidate_K < best_K ||
              (candidate_K == best_K && candidate_Q < best_Q)) {
            best_Q = candidate_Q;
            best_K = candidate_K;
            best_s = s;
          }

          if (write != k) {
            index[write] = s;
            tests[write] = std::move(tests[k]);
          }
          last_updates[write] = last_up;
          ++write;
        }
      }

      index.resize(write);
      tests.erase(tests.begin() + write, tests.end());
      last_updates.resize(write);
    } else {
      for (size_t k = 0; k < m; ++k) {
        const size_t s = index[k];

        for (size_t u = last_updates[k] + 1; u <= t; ++u) {
          tests[k].update(data[u - 1]);
        }
        last_updates[k] = t;

        if (tests[k].statistic() < gamma) {
          const size_t candidate_K = K[s] + 1;
          const double candidate_Q = Q[s] + segment_cost(S1, S2, s, t);

          if (candidate_K < best_K ||
              (candidate_K == best_K && candidate_Q < best_Q)) {
            best_Q = candidate_Q;
            best_K = candidate_K;
            best_s = s;
          }
        }
      }
    }

    Q[t] = best_Q;
    K[t] = best_K;
    previous[t] = best_s;

    index.push_back(t);
    tests.emplace_back();
    last_updates.push_back(t);
  }

  std::vector<size_t> changepoints;
  changepoints.reserve(n / 2 + 1);
  for (size_t i = n; i > 0; i = previous[i]) {
    changepoints.push_back(i);
  }
  std::reverse(changepoints.begin(), changepoints.end());
  std::reverse(index.begin(), index.end());

  return List::create(
    _["changepoints"] = changepoints,
    _["lastIndexSet"] = index,
    _["nb"] = nb,
    _["costQ"] = NULL,
    _["R"] = build_R_matrix(Q, K, previous)
  );
}

} // namespace

// [[Rcpp::export]]
List SVP(std::vector<double> data,
         double gamma,
         std::string test,
         bool prune_if_unvalid = true)
{
  if (test == "gaussian_mean") {
    return svp_impl<GaussianMean>(data, gamma, prune_if_unvalid);
  }
  if (test == "gamma_rate") {
    return svp_impl<GammaRate>(data, gamma, prune_if_unvalid);
  }
  if (test == "gaussian_variance") {
    return svp_impl<GaussianVariance>(data, gamma, prune_if_unvalid);
  }

  stop("Unknown test type");
}
