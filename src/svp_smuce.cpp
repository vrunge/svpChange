#include <Rcpp.h>
using namespace Rcpp;

//' C++ SVP with SMUCE validity and constrained Gaussian cost
//' @param y Numeric observations.
//' @param q SMUCE threshold.
//' @param sigma2 Known Gaussian variance.
//' @return Integer segment-end indices.
//' @details
//' For independent Gaussian observations with known variance, this function
//' returns the segment endpoints of one SMUCE-optimal partition. Internal
//' endpoints are estimated changepoints.
//' @references
//' Frick, K., Munk, A., and Sieling, H. (2014). Multiscale Change-Point
//' Inference. *Journal of the Royal Statistical Society: Series B*, 76(3),
//' 495--580. doi:10.1111/rssb.12047.
//' @examples
//' set.seed(1)
//' data <- ts_generator(
//'   chpts = c(20, 40), parameters = c(0, 2),
//'   sd_noise = 1, type = "gauss"
//' )
//' svp_smuce_cpp(data, q = 1.5, sigma2 = 1)
// [[Rcpp::export]]
IntegerVector svp_smuce_cpp(NumericVector y, double q, double sigma2 = 1.0) {
  int n = y.size();
  if (n < 1 || !R_finite(q) || q < 0 || !R_finite(sigma2) || sigma2 <= 0)
    stop("invalid arguments");
  std::vector<double> cs(n + 1, 0.0), cs2(n + 1, 0.0);
  for (int i = 0; i < n; ++i) {
    if (!R_finite(y[i])) stop("y must contain only finite values");
  }
  // Center before prefix sums to avoid cancellation for translated data.
  const double offset = y[0];
  for (int i = 0; i < n; ++i) {
    const double centered = y[i] - offset;
    cs[i+1] = cs[i] + centered;
    cs2[i+1] = cs2[i] + centered * centered;
  }

  // The SMUCE radius depends only on the interval length.
  std::vector<double> radius(n + 1, 0.0);
  for (int len = 1; len <= n; ++len) {
    double length = static_cast<double>(len);
    radius[len] = std::sqrt(sigma2 / length) *
      (q + std::sqrt(2.0 * (1.0 + std::log(n / length))));
  }

  std::vector<int> K(n+1, n+1), prev(n+1, -1);
  std::vector<double> C(n+1, R_PosInf);
  K[0] = 0; C[0] = 0.0;

  // lo_previous[s] and hi_previous[s] describe the admissible mean interval
  // for segment [s, t - 1]. The current column is computed for decreasing s
  // from the two neighbouring segments and the local interval [s, t].
  std::vector<double> lo_previous(n, R_NegInf);
  std::vector<double> hi_previous(n, R_PosInf);
  std::vector<double> lo_current(n, R_NegInf);
  std::vector<double> hi_current(n, R_PosInf);

  for (int t = 0; t < n; ++t) {
    int first_valid_start = 0;

    for (int s = t; s >= 0; --s) {
      int len = t - s + 1;
      double segment_sum = cs[t+1] - cs[s];
      double mean_segment = segment_sum / len;
      double local_lo = mean_segment - radius[len];
      double local_hi = mean_segment + radius[len];

      if (s == t) {
        lo_current[s] = local_lo;
        hi_current[s] = local_hi;
      } else {
        lo_current[s] = std::max(
          local_lo,
          std::max(lo_current[s+1], lo_previous[s])
        );
        hi_current[s] = std::min(
          local_hi,
          std::min(hi_current[s+1], hi_previous[s])
        );
      }

      // If [s, t] is invalid, every [s', t] with s' < s is invalid too.
      if (lo_current[s] > hi_current[s]) {
        first_valid_start = s + 1;
        break;
      }
    }

    // Process starts in increasing order to preserve the previous tie rule.
    for (int s = first_valid_start; s <= t; ++s) {
      if (K[s] > n || !R_finite(C[s])) continue;

      int len = t - s + 1;
      double segment_sum = cs[t+1] - cs[s];
      double mean_segment = segment_sum / len;
      double theta = std::min(
        std::max(mean_segment, lo_current[s]),
        hi_current[s]
      );
      double rss = (cs2[t+1] - cs2[s]) -
        2.0 * theta * segment_sum + len * theta * theta;
      rss = std::max(0.0, rss);
      double cost = C[s] + rss / sigma2;
      int nk = K[s] + 1;
      if (nk < K[t+1] || (nk == K[t+1] && cost < C[t+1])) {
        K[t+1] = nk;
        C[t+1] = cost;
        prev[t+1] = s;
      }
    }

    std::swap(lo_previous, lo_current);
    std::swap(hi_previous, hi_current);
  }
  std::vector<int> ends; int t = n;
  while (t > 0) {
    if (prev[t] < 0) stop("no valid partition");
    ends.push_back(t);
    t = prev[t];
  }
  std::reverse(ends.begin(), ends.end());
  return IntegerVector(ends.begin(), ends.end());
}
