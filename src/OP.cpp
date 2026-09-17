#include <Rcpp.h>
#include <vector>
#include <limits>
#include <algorithm>

using namespace Rcpp;

//' Optimal Partitioning Algorithm
//'
//' @title Optimal Partitioning Algorithm
//' @description Finds the least-squares segmentation minimizing the sum of
//' within-segment squared errors plus `penalty` times the number of estimated
//' change points.
//'
//' @param data A numeric vector representing the data to segment.
//' @param penalty Numeric penalty applied to each estimated change point.
//'
//' @details A candidate boundary `s` and endpoint `t` represent the R segment
//' `data[(s + 1):t]`. Setting the initial cost to `-penalty` makes the total
//' penalty equal to `penalty * (K - 1)` for a partition with `K` segments.
//'
//' @return A list with the following components:
//' \describe{
//'   \item{changepoints}{Increasing, one-based, inclusive segment endpoints,
//'     including `length(data)`.}
//'   \item{lastIndexSet}{Always `NULL`; OP does not prune candidates.}
//'   \item{nb}{Always `NULL`; OP does not prune candidates.}
//'   \item{costQ}{Numeric vector of length `length(data)`. Element `t` is the
//'     minimum penalized cost for `data[1:t]`.}
//' }
//'
//' @examples
//' set.seed(1)
//' data <- ts_generator(
//'   chpts = c(40, 80, 120), parameters = c(0, 2, -1),
//'   sd_noise = 1, type = "gauss"
//' )
//' penalty <- 2 * log(length(data))
//' OPres <- OP(data, penalty)
//' OPres$changepoints
//'
//' @seealso [PELT()] for the pruned version of the same objective, [SN()] for
//'   fixed numbers of segments, and [SVP()] for validity-constrained
//'   partitioning.
//' @export
// [[Rcpp::export]]
List OP(std::vector<double> data, double penalty)
{
  size_t n = data.size();

  // Initialize the costs and the changepoints
  std::vector<double> Q(n + 1, std::numeric_limits<double>::infinity());
  Q[0] = -penalty;
  std::vector<size_t> lastChange(n + 1, 0);

  //
  // PREPROCESSING
  //
  // Cumulative sum for optimized calculations
  std::vector<double> S1(n + 1, 0);
  std::vector<double> S2(n + 1, 0);

  for (size_t i = 0; i < n; i++)
  {
    S1[i + 1] = S1[i] + data[i];
    S2[i + 1] = S2[i] + data[i] * data[i];
  }

  //
  // MAIN LOOP
  //
  double tempQ;

  // Cost calculation for each sub-segment
  for (size_t t = 1; t <= n; t++)
  {
    for (size_t s = 0; s < t; s++)
    {
      // Segment cost [s+1, t]
      // Total cost with beta penalization
      tempQ = Q[s] + (S2[t] - S2[s]) - (S1[t] - S1[s]) * (S1[t] - S1[s]) / (t - s) + penalty;

      // Minimization
      if (tempQ < Q[t])
      {
        Q[t] = tempQ;
        lastChange[t] = s;
      }
    }
  }

  //
  // BACKTRACKING
  //
  // Change points reconstruction
  std::vector<int> changepoints;
  size_t i = n;
  while (lastChange[i] > 0)
  {
    changepoints.push_back(static_cast<int>(lastChange[i]));
    i = lastChange[i];
  }
  std::reverse(changepoints.begin(), changepoints.end());
  changepoints.push_back(static_cast<int>(n));

  //
  // Return
  //
  return List::create(
    Named("changepoints") = changepoints,
    Named("lastIndexSet") = R_NilValue,
    Named("nb")           = R_NilValue,
    Named("costQ")        = std::vector<double>(Q.begin() + 1, Q.end()));
}
