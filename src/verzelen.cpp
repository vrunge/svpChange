#include "svp_verzelen.h"

// Internal exact statistic, with an exhaustive switch for regression checks.
// [[Rcpp::export(name = ".verzelen_statistic_cpp")]]
double verzelen_statistic_cpp(std::vector<double> data, double sigma2 = 1.0,
                             int N = 0, bool accelerated = true)
{
  if (data.empty() || !std::all_of(data.begin(), data.end(),
      [](double x) { return std::isfinite(x); }))
    Rcpp::stop("'data' must contain finite observations and be nonempty.");
  if (!std::isfinite(sigma2) || sigma2 <= 0.0)
    Rcpp::stop("sigma2 must be finite and positive");
  if (N == 0) N = static_cast<int>(data.size());
  if (N < static_cast<int>(data.size()))
    Rcpp::stop("N must be at least the data length");
  return svp_detail::verzelen_statistic(data, sigma2, N, accelerated);
}
