#include <Rcpp.h>
#include <vector>
#include <limits>
#include <algorithm>

using namespace Rcpp;


//' Segment Neighborhood
//'
//' @title Segment Neighborhood
//'
//' @description This function implements the SN algorithm of a given vector `data` with a given maximum number of changes
//' It finds the optimal segmentation for all K between 1 and Kmax that minimizes the global cost using dynamic programming.
//'
//' @param data A numeric vector representing the data to segment.
//' @param Kmax An integer value representing the maximal number of segment
//'
//' @return A list with the following elements:
//' \itemize{
//'   \item \code{changepoints}: the last index of each segment for each number of segment K between 1 and Kmax
//'   \item \code{nb}: NULL
//'   \item \code{lastIndexSet}: NULL
//'   \item \code{costQ}: a matrix saving the optimal cost at each time step for each number of segments (in k-th column for k segments)
//' }
//'
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
    for (size_t t = 1; t < n + 1; t++)
    {
      for (size_t s = 0; s < t; s++)
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
  IntegerVector cp(1);   // k segments
  cp[0] = n;
  changepoints[0] = cp;

  for (int k = 1; k < Kmax; k++)
  {
    IntegerVector cp(k+1);   // k segments
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
    Named("lastIndexSet") = NULL,
    Named("nb") = NULL,
    Named("costQ") = costQ);
}










