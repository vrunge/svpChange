#include <Rcpp.h>
#include <vector>

using namespace Rcpp;

//' Smallest Valid Partitioning with Validation and Pruning using Rcpp
//'
//' @title Smallest Valid Partitioning with Validation and Pruning
//' @description This function implements a dynamic programming approach to segment a univariate signal into the smallest number of valid segments, according to a user-defined validation function. Each segment must pass a validity test (e.g., based on variance, range, etc.). The algorithm minimizes a quadratic cost subject to this constraint.
//'
//' @param data A numeric vector representing the univariate signal to be segmented.
//' @param gamma A numeric value used as a threshold in the validation function and as a penalty for each segment.
//' @param test A function of the form `function(data, gamma)` returning TRUE if the segment is valid. Default is `valid_OP`.
//' @param all_full_validity Logical. If TRUE (default), the algorithm applies *segment-wise validation*:
//' at each time step, it tests whether the candidate segment \code{data[(s+1):t]} is valid using the
//' user-defined function \code{test}. If the segment fails the test, the candidate \code{s} is removed
//' (pruned) from the set of possible changepoints. This accelerates computation by avoiding invalid
//' segment extensions. If FALSE, the algorithm skips this validation and considers all candidate
//' segments without checking their validity (which can be faster but may return invalid segments).
//'
//' @return A list with the following components :
//' \describe{
//'   \item{changepoints}{Integer vector indicating the ending index of each segment (i.e., positions of changepoints).}
//'   \item{nb}{Integer vector of length \code{length(data)}. At each position \code{t}, it records the number of candidates tested.}
//'   \item{costQ}{Numeric vector of length \code{length(data)}. Quadratic cost value at each time step.}
//'   \item{R}{A matrix of dimension \code{(length(data)+1) x 3} containing, for each time step :
//'     \describe{
//'       \item{Q}{cumulative cost}
//'       \item{K}{number of segments}
//'       \item{s}{previous changepoint}
//'     }
//'   }
//' }
//'
//' @export
// [[Rcpp::export]]
List svp0(std::vector<double> data,
          double gamma,
          Function test,
          bool all_full_validity = true)
{
  size_t n = data.size();

   // Initialization of the R matrix
   NumericMatrix R(n + 1, 3); // Q, K, s
   std::fill(R.begin(), R.end(), R_PosInf);
   R(0, 0) = 0.0; // Default value
   R(0, 1) = 0.0; // Default value
   R(0, 2) = 0.0; // Default value

   IntegerVector nb(n); // nb de candidats examinés à chaque t
   NumericVector costQ(n); // cost vector

   // Cumulative sum for optimized calculations
   NumericVector cs_y(n + 1, 0.0);
   NumericVector cs_y2(n + 1, 0.0);
   for (int i = 0; i < n; ++i)
   {
     cs_y[i + 1] = cs_y[i] + data[i];
     cs_y2[i + 1] = cs_y2[i] + data[i] * data[i];
   }

   std::vector<int> INDEX = {0};


   bool valid;

   for (int t = 1; t <= n; t++)
   {
     double best_Q = R_PosInf;
     int best_K = INT_MAX;
     int best_s = -1;

     std::vector<int> new_INDEX;
     for (size_t k = 0; k < INDEX.size(); ++k)
     {
       int s = INDEX[k];
       if (s >= t) continue;

      valid = true;
      if (all_full_validity)
      {
        std::vector<double> seg(data.begin() + s, data.begin() + t);
        valid = as<bool>(test(seg, gamma));
      }

      if (valid == true)
      {
        double candidate_Q = R(s, 0) + (cs_y2[t] - cs_y2[s]) - ((cs_y[t] - cs_y[s]) * (cs_y[t] - cs_y[s])) / (t - s);
        int candidate_K = R(s, 1) + 1;

        if (candidate_K < best_K || (candidate_K == best_K && candidate_Q < best_Q))
        {
          best_Q = candidate_Q;
          best_K = candidate_K;
          best_s = s;
          costQ[t - 1] = best_Q;
        }
        new_INDEX.push_back(s);
      }
    }

    //
    // PRUNING PELT
    //
    std::vector<int> pruned_INDEX;
    for (size_t k = 0; k < new_INDEX.size(); ++k)
    {
      int s = new_INDEX[k];
      if (s == t) continue;

      double candidate_Q = R(s, 0) + (cs_y2[t] - cs_y2[s]) - ((cs_y[t] - cs_y[s]) * (cs_y[t] - cs_y[s])) / (t - s);
      int candidate_K = R(s, 1);

      if (!(candidate_Q > best_Q && candidate_K == best_K))
      {
        pruned_INDEX.push_back(s);
      }
    }
    pruned_INDEX.push_back(t);
    INDEX = pruned_INDEX;
    nb[t - 1] = INDEX.size();

    R(t, 0) = best_Q;
    R(t, 1) = best_K;
    R(t, 2) = best_s;
  }

  // Reconstruct changepoints
  IntegerVector changepoints;
  int t = n;
  while (t > 0)
  {
    changepoints.push_front(t);
    t = R(t, 2);
  }

  return List::create(
    _["changepoints"] = changepoints,
    _["nb"] = nb,
    _["costQ"] = costQ,
    _["R"] = R
  );
}
