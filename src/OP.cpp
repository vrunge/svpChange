#include <Rcpp.h>
#include <vector>
#include <limits>
#include <algorithm>

using namespace Rcpp;

//' Optimal Partitioning Algorithm
//'
//' @title Optimal Partitioning Algorithm
//' @description This function implements the OP algorithm of a given vector `data` with a given penalty term.
//' It finds the optimal change points that minimize a penalized cost using dynamic programming.
//'
//' @param data A numeric vector representing the data to segment.
//' @param penalty A double value representing the penalty term for adding a new segment.
//'
//' @return A list with (1) the `\code{changepoints}` elements, (2) a vector `\code{nb}` saving the number of non-pruned elements at each iteration, (3) a vector `\code{lastIndexSet}` containing the non-pruned indices at the end of the algo and (4) a vector `\code{costQ}` saving the optimal cost at each time step. The elements `\code{nb}` and `\code{lastIndexSet}` are set to NULL. They are included for consistency with other algorithms that prune the number of indices to be tracked over time.
//'
//' @examples
//' n <- 1000
//' data <- rep(c(0, 1, -0.5, 0), each = n) + rnorm(4 * n)
//' penalty <- 2 * log(length(data))
//' OPres <- OP(data, penalty)
//' OPres$changepoints
//'
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
  std::vector<size_t> changepoints;
  size_t i = n;
  while (lastChange[i] > 0)
  {
    changepoints.push_back(lastChange[i]);
    i = lastChange[i];
  }
  std::reverse(changepoints.begin(), changepoints.end());
  changepoints.push_back(n);

  //
  // Return
  //
  return List::create(
    Named("changepoints") = changepoints,
    Named("lastIndexSet") = NULL,
    Named("nb")           = NULL,
    Named("costQ")        = std::vector<double>(Q.begin() + 1, Q.end()));
}

