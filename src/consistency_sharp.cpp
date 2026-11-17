#include "svp_consistency_sharp.h"

// Internal exact statistic, for calibration and regression checks.
// [[Rcpp::export(name = ".consistency_sharp_statistic_cpp")]]
double consistency_sharp_statistic_cpp(std::vector<double> data, double sigma2 = 1.0,
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
  return svp_detail::consistency_sharp_statistic(data, sigma2, N, accelerated);
}
