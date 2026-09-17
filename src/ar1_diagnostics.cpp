#include "svp_ar1_helpers.h"
#include "svp_tests_ar1.h"
#include <Rcpp.h>
#include <cmath>
#include <vector>

using namespace Rcpp;

// [[Rcpp::export]]
double AR1_rho(std::vector<double> data)
{
  return svp_detail::robust_ar1_rho(data);
}

// [[Rcpp::export]]
List AR1_single_change(std::vector<double> data,
                       double gamma,
                       double rho = NA_REAL,
                       double sigma2 = 1.0,
                       bool profile_sigma = false)
{
  const double rho_used = NumericVector::is_na(rho) ?
    svp_detail::robust_ar1_rho(data) : rho;
  if (!std::isfinite(rho_used) || std::fabs(rho_used) >= 1.0) {
    stop("rho must be finite and strictly between -1 and 1");
  }
  if (!profile_sigma &&
      (!std::isfinite(sigma2) || sigma2 <= 0.0)) {
    stop("sigma2 must be finite and positive");
  }

  AR1ExactMeanChange test(rho_used, sigma2, profile_sigma);
  for (double value : data) test.update(value);
  const double statistic = test.statistic();
  const int changepoint = test.changepoint();

  return List::create(
    _["rho"] = rho_used,
    _["sigma2"] = sigma2,
    _["profile_sigma"] = profile_sigma,
    _["statistic"] = statistic,
    _["rss0"] = test.rss0(),
    _["rss1"] = test.rss1(),
    _["changepoint"] = changepoint > 0 ? changepoint : NA_INTEGER,
    _["valid"] = statistic < gamma
  );
}
