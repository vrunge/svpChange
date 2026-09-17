#include <Rcpp.h>
#include <vector>
#include <limits>
#include <algorithm>

using namespace Rcpp;

//' Segment Neighborhood
//'
//' @title Segment Neighborhood
//'
//' @description Finds the least-squares segmentation for every fixed number of
//' segments from 1 through `Kmax`, using dynamic programming.
//'
//' @param data A numeric vector representing the data to segment.
//' @param Kmax Maximum number of segments. It should be an integer between 1
//'   and `length(data)`.
//'
//' @details A candidate boundary `s` and endpoint `t` represent the R segment
//' `data[(s + 1):t]`. The segment cost is its sum of squared errors around its
//' sample mean.
//'
//' @return A list with the following components:
//' \describe{
//'   \item{changepoints}{List of length `Kmax`. Element `k` contains the `k`
//'     increasing, one-based, inclusive segment endpoints of the optimal
//'     `k`-segment partition, including `length(data)`.}
//'   \item{lastIndexSet}{Always `NULL`; SN does not prune candidates.}
//'   \item{nb}{Always `NULL`; SN does not prune candidates.}
//'   \item{costQ}{Numeric matrix with `length(data) + 1` rows and `Kmax`
//'     columns. Entry `[t + 1, k]` is the minimum cost for `data[1:t]` with
//'     exactly `k` segments; infeasible entries are `Inf`.}
//' }
//'
//' @examples
//' set.seed(1)
//' data <- ts_generator(
//'   chpts = c(20, 40), parameters = c(0, 4),
//'   sd_noise = 0.5, type = "gauss"
//' )
//' fit <- SN(data, Kmax = 2)
//' fit$changepoints[[2]]
//' fit$costQ[nrow(fit$costQ), ]
//'
//' @seealso [OP()] and [PELT()] for penalized segmentation, and [SVP()] for
//'   validity-constrained partitioning.
//' @export
// [[Rcpp::export]]
List SN(std::vector<double> data, int Kmax)
{
  size_t n = data.size();

  NumericMatrix costQ(n + 1, Kmax);      // C[t,k]
  IntegerMatrix lastChange(n + 1, Kmax); // argmin s for (t,k), store s (0..t-1)

  // Initialize with +inf
  for (int t = 0; t < n + 1; t++)
    for (int k = 0; k < Kmax; k++)
    {
      costQ(t, k) = std::numeric_limits<double>::infinity();
      lastChange(t, k) = 0;
    }


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

  // base case K = 1
  for (int t = 1; t < n + 1; t++)
  {
    costQ(t, 0) = S2[t] - S1[t] * S1[t] / t;
    lastChange(t, 0) = 0;
  }

  //
  // MAIN LOOP
  //
  double tempQ;

  // Cost calculation for each sub-segment
  for (int k = 1; k < Kmax; k++)
  {
    const size_t minimum_start = static_cast<size_t>(k);
    for (size_t t = minimum_start + 1; t < n + 1; t++)
    {
      for (size_t s = minimum_start; s < t; s++)
      {
        tempQ = costQ(s,k-1) + (S2[t] - S2[s]) - (S1[t] - S1[s]) * (S1[t] - S1[s]) / (t - s);
        if (tempQ < costQ(t,k))
        {
          costQ(t,k) = tempQ;
          lastChange(t,k) = s;
        }
      }
    }
  }

  //
  // BACKTRACKING
  //
  // Change points reconstruction
  List changepoints(Kmax);
  IntegerVector cp(1);   // one segment
  cp[0] = n;
  changepoints[0] = cp;

  for (int k = 1; k < Kmax; k++)
  {
    IntegerVector cp(k+1);   // k + 1 segments
    cp[k] = n;
    int t = n;
    int curr_k = k;

    for (int i = (k - 1); i >= 0; --i)
    {
      const int s = lastChange(t, curr_k);
      cp[i] = s;          // changepoint at position s (end of previous segment)
      t = s;              // continue reconstructing up to s
      curr_k = curr_k - 1;
    }
    changepoints[k] = cp;
  }

  return List::create(
    Named("changepoints") = changepoints,
    Named("lastIndexSet") = R_NilValue,
    Named("nb") = R_NilValue,
    Named("costQ") = costQ);
}







