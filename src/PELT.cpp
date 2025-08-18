#include <Rcpp.h>
#include <vector>
#include <limits>
#include <algorithm>

using namespace Rcpp;

//' Optimal Partitioning algorithm using PELT
//'
//' @title Optimal Partitioning using PELT
//' @description This function implements the OP algorithm using PELT of a given vector `data` with a given penalty term.
//' It finds the optimal changepoints that minimize the cost function using dynamic programming.
//'
//' @param data A numeric vector representing the data to segment.
//' @param penalty A double value representing the penalty term for adding a new segment.
//'
//' @return A list with
//' (1) the changepoint elements (each last index of each segment in \code{changepoints}),
//' (2) a vector \code{nb} saving the number of non-pruned elements at each iteration,
//' (3) a vector \code{lastIndexSet} containing the non-pruned indices at the end of the algorithm,
//' (4) a vector \code{costQ} saving the optimal cost at each time step.
//'
//' @examples
//' n <- 1000
//' data <- rep(c(0, 5, 2.5, 7), each = n) + rnorm(4 * n)
//' penalty <- 2 * log(length(data))
//' resPELT <- PELT(data, penalty)
//'
//' @export
// [[Rcpp::export]]
List PELT(std::vector<double> data, double penalty)
{
  size_t n = data.size();

  // Initialize the costs and the changepoints
  std::vector<double> Q(n + 1, std::numeric_limits<double>::infinity());
  Q[0] = -penalty;
  std::vector<size_t> last_cp(n + 1, 0);
  std::vector<size_t> P(1, 0);
  std::vector<size_t> length_P(n);

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
  double best_cost = std::numeric_limits<double>::infinity();
  double sum_x;
  double sum_x2;
  double gaussian_cost;
  size_t t1;
  size_t arg_min;
  size_t s;

  for (size_t t = 0; t < n; t++)
  {
    t1 = t + 1;
    std::vector<double> costs(P.size(), std::numeric_limits<double>::infinity());
    best_cost = std::numeric_limits<double>::infinity();
    arg_min = 1;

    for (size_t i = 0; i < P.size(); i++)
    {
      s = P[i];
      sum_x = S1[t1] - S1[s];
      sum_x2 = S2[t1] - S2[s];
      gaussian_cost = sum_x2 - (sum_x * sum_x) / (t1 - s);
      costs[i] = Q[s] + gaussian_cost + penalty;
      if (costs[i] < best_cost)
      {
        best_cost = costs[i];
        arg_min = s;
      }
    }
    Q[t1] = best_cost;
    last_cp[t1] = arg_min;
    length_P[t] = P.size();

    //
    // PELT PRUNING
    //
    std::vector<size_t> newP;
    for (size_t i = 0; i < P.size(); i++)
    {
      size_t s = P[i];
      if (costs[i] <= Q[t1] + penalty) newP.push_back(s);
    }
    newP.push_back(t1);
    P = newP;
  }

  //
  // BACKTRACKING
  //
  // Change points reconstruction
  std::vector<int> changepoints;
  size_t i = n;
  while (last_cp[i] > 0)
  {
    changepoints.push_back(last_cp[i]);
    i = last_cp[i];
  }
  std::reverse(changepoints.begin(), changepoints.end());
  changepoints.push_back(n);

  std::reverse(P.begin(), P.end());

  return List::create(
    Named("changepoints") = wrap(changepoints),
    Named("lastIndexSet") = wrap(P),
    Named("nb") = wrap(length_P),
    Named("costQ") = wrap(std::vector<double>(Q.begin() + 1, Q.end())));
}
