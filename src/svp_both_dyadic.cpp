#include <Rcpp.h>

#include "gaussian_focus_state.h"
#include "svp_internal.h"
#include "svp_test_base.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

namespace {

std::size_t dyadic_band_count(std::size_t n)
{
  std::size_t count = 0;
  std::size_t upper = n / 2;
  while (upper >= 2) {
    ++count;
    upper /= 2;
  }
  return count + 1;
}

std::size_t dyadic_band(std::size_t n, std::size_t length)
{
  const std::size_t count = dyadic_band_count(n);
  std::size_t band = 0;
  std::size_t lower = n / 2;
  while (band + 1 < count && length <= lower) {
    ++band;
    lower /= 2;
  }
  return band;
}

class GaussianMeanDyadic : public TestBase
{
public:
  GaussianMeanDyadic(const std::vector<double>& gamma, std::size_t n)
    : gamma_(gamma), n_(n), length_(0) {}

  void update(double y) override
  {
    info_.update(y);
    ++length_;
  }

  double statistic() const override
  {
    return info_.max_cost();
  }

  bool passes(double) const override
  {
    const double threshold = gamma_[dyadic_band(n_, length_)];
    return !info_.exceeds_cost(threshold);
  }

private:
  std::vector<double> gamma_;
  std::size_t n_;
  std::size_t length_;
  changepoint::GaussianFocusState info_;
};

} // namespace

//' Smallest Valid Partitioning with dyadic Gaussian thresholds and both pruning
//'
//' The threshold used for a candidate segment depends on its length. For a
//' series of length `n`, band `j` covers lengths in
//' `(n / 2^(j + 1), n / 2^j]`, with `j = 0` denoting the longest segments.
//' Single-observation segments are always valid. Candidate removal follows
//' the same both-direction pruning rules as [SVP()].
//'
//' @param data Finite numeric observations.
//' @param gamma Finite positive numeric vector with one threshold per dyadic
//'   band, ordered from longest to shortest segment lengths.
//' @return The same partition and dynamic-programming diagnostics returned by
//'   [SVP()].
//' @export
// [[Rcpp::export]]
Rcpp::List SVP_both_dyadic(std::vector<double> data,
                           Rcpp::NumericVector gamma)
{
  if (data.empty()) {
    Rcpp::stop("'data' must contain at least one observation.");
  }
  if (!std::all_of(data.begin(), data.end(),
                   [](double value) { return std::isfinite(value); })) {
    Rcpp::stop("'data' must contain only finite observations.");
  }
  if (gamma.size() != static_cast<R_xlen_t>(dyadic_band_count(data.size()))) {
    Rcpp::stop("'gamma' must contain one value per dyadic length band.");
  }
  if (!std::all_of(gamma.begin(), gamma.end(),
                   [](double value) { return std::isfinite(value) && value > 0.0; })) {
    Rcpp::stop("'gamma' values must be finite and positive.");
  }

  const std::vector<double> thresholds(gamma.begin(), gamma.end());
  const svp_detail::SegmentCost cost =
    svp_detail::SegmentCost::gaussian(data);

  // svp_impl() supplies the established both-pruning and result semantics.
  // The dyadic test selects its threshold from its incremental segment length.
  return svp_detail::svp_impl<GaussianMeanDyadic>(
    data, thresholds.front(), "both", cost, thresholds, data.size()
  );
}
