#include <Rcpp.h>
#include <vector>
#include <limits>
#include <algorithm>

using namespace Rcpp;

//' Segment neighborhood using PELT
//'
//' @title Segment neighborhood using PELT
//'
//' @description This function implements the SN algorithm using PELT of a given vector `data` with a given maximum number of changes
//' It finds the optimal segmentation for all K between 1 and Kmax that minimize the global cost using dynamic programming.
//'
//' @param data A numeric vector representing the data to segment.
//' @param Kmax An integer value representing the
//'
//' @return A list with the following elements:
//' \itemize{
//'   \item \code{changepoints}: the last index of each segment,
//'   \item \code{nb}: a vector saving the number of non-pruned elements at each iteration,
//'   \item \code{lastIndexSet}: a vector containing the non-pruned indices at the end of the algorithm,
//'   \item \code{costQ}: a vector saving the optimal cost at each time step.
//' }
//'
//' @export
// [[Rcpp::export]]
List SN(std::vector<double> data, int Kmax)
{
  size_t n = data.size();


  return List::create(
    Named("changepoints") = NULL,
    Named("lastIndexSet") = NULL,
    Named("nb") = NULL,
    Named("costQ") = NULL);
}










