#include "svp_ar1_helpers.h"
#include "svp_internal.h"
#include "svp_verzelen.h"
#include "svp_consistency_sharp.h"
#include "svp_tests_empirical_both.h"
#include "validity_tests.h"
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>

using namespace Rcpp;

// -----------------------------------------------------------------------------
// Test dispatch
// -----------------------------------------------------------------------------

namespace {

// Names accepted by the public SVP() dispatcher.
enum class TestKind {
  GaussianMean,
  GammaRate,
  GaussianVariance,
  Quantile,
  QuantileExact,
  VarCost,
  WilcoxonCost,
  MedianMoodCost,
  AR1,
  AR1Profile,
  AR1Focus
};

// Within-segment cost used to rank partitions with the same segment count.
enum class CostKind {
  Auto,
  Gaussian,
  AR1
};

// Convert the user-facing cost name to the corresponding internal enum value.
CostKind parse_cost_kind(const std::string& cost)
{
  if (cost == "auto") return CostKind::Auto;
  if (cost == "gaussian") return CostKind::Gaussian;
  if (cost == "ar1") return CostKind::AR1;

  stop("Unknown cost type: '%s'.", cost.c_str());
  return CostKind::Auto;
}

// Convert the user-facing test name to the corresponding internal enum value.
TestKind parse_test_kind(const std::string& test)
{
  if (test == "gaussian_mean") return TestKind::GaussianMean;
  if (test == "gamma_rate") return TestKind::GammaRate;
  if (test == "gaussian_variance") return TestKind::GaussianVariance;
  if (test == "quantile") return TestKind::Quantile;
  if (test == "quantileExact") return TestKind::QuantileExact;
  if (test == "varCost") return TestKind::VarCost;
  if (test == "WilcoxonCost") return TestKind::WilcoxonCost;
  if (test == "MedianMoodCost") return TestKind::MedianMoodCost;
  if (test == "AR1") return TestKind::AR1;
  if (test == "AR1Profile") return TestKind::AR1Profile;
  if (test == "AR1Focus") return TestKind::AR1Focus;

  stop("Unknown test type: '%s'.", test.c_str());
  return TestKind::GaussianMean;
}

} // namespace

// -----------------------------------------------------------------------------
// Public SVP API
// -----------------------------------------------------------------------------

//' Smallest Valid Partitioning with Incremental Validity Tests
//'
//' Segments a univariate series into the smallest number of segments that pass
//' a selected validity test. Among partitions with the same number of segments,
//' the function chooses the one with the smallest within-segment cost, selected
//' with `cost`. The validity statistics are maintained incrementally in C++, so
//' this interface is considerably faster than supplying an R validity function
//' to [svp0()].
//'
//' A candidate segment is evaluated only at its current endpoint. The
//' incremental state still receives every intervening observation, but a
//' previously evaluated endpoint is not tested again. Whether an invalid
//' endpoint can be reconsidered later depends on `subtests`: with
//' `subtests = "none"`, the candidate is retained and can recover; with
//' `subtests = "right"` or `subtests = "both"`, it is removed immediately.
//' With `subtests = "both"`, an invalid boundary also removes that boundary
//' and every older boundary under the inclusive left-pruning rule. The
//' available values of `test` are:
//'
//' * `"gaussian_mean"`: Gaussian FOCUS likelihood-ratio test for a change in
//'   mean. Use this for independent Gaussian observations with constant
//'   variance.
//' * `"gamma_rate"`: likelihood-ratio test for a change in the rate of
//'   positive Gamma observations (unit shape in the current implementation).
//' * `"gaussian_variance"`: Gamma-rate test applied to squared observations,
//'   for changes in Gaussian variance around a known zero mean.
//' * `"AR1"`: conditional fixed-`rho` Gaussian likelihood-ratio scan for a
//'   change in the marginal mean of an AR(1) series, conditional on the first
//'   observation. In the AR(1) model, an innovation is
//'   the new random shock after accounting for the previous observation and the
//'   AR(1) mean structure; `sigma2` is the variance of this shock, not the
//'   marginal variance of the observations.
//' * `"AR1Profile"`: the same conditional AR(1) scan, but profiling out the
//'   innovation variance. The variance is estimated separately under the no-change and
//'   change models from their residual sums of squares, which is useful when
//'   the innovation scale is unknown.
//' * `"AR1Focus"`: faster approximate AR(1) test that applies Gaussian FOCUS
//'   to the innovations `x[t] - rho * x[t - 1]`. The exact `"AR1"` test is the
//'   preferred choice when boundary accuracy matters.
//' * `"quantile"`: approximate quantile-range test using two streaming
//'   P-squared quantile estimators.
//' * `"quantileExact"`: exact empirical quantile-range test using an
//'   order-statistics tree.
//' * `"varCost"`: online within-segment variance statistic.
//' * `"WilcoxonCost"`: maximum Wilcoxon statistic over possible splits.
//' * `"MedianMoodCost"`: maximum Median-Mood statistic over possible splits.
//' For Gaussian means, `subtests = "both_msl"`, `"both_ms_refined"`, and
//' `"both_focus_lil"` apply the `"both"` pruning rule to multiscale tests.
//'
//' ## Per-candidate update complexity
//'
//' Let `m` be the current length of one candidate segment and let `a_m` be
//' the number of active functional-pruning pieces, where `a_m <= m`. The
//' state-update column is the cost of incorporating one new observation.
//' The endpoint-check column includes the work needed to decide validity.
//' These are per-candidate costs; total SVP time also depends on how many
//' candidates are examined under the selected `subtests` mode.
//'
//' | Test | State update | Update and endpoint validity check |
//' |---|---:|---:|
//' | `"gaussian_mean"` | amortized `O(1)`, worst `O(m)` | expected `O(log m)`, worst `O(m)` |
//' | `"gaussian_mean"`, `"both_focus_lil"` | amortized `O(1)`, worst `O(m)` | exact FOCUS maximum, expected `O(log m)`, worst `O(m)` |
//' | `"gaussian_mean"`, `"both_msl"` or `"both_ms_refined"` | amortized `O(1)`, worst `O(m)` | FOCUS bound, exact `O(m)` scan if needed |
//' | `"gamma_rate"` | `O(a_m)`, worst `O(m)` | `O(a_m)`, worst `O(m)` |
//' | `"gaussian_variance"` | `O(a_m)`, worst `O(m)` | `O(a_m)`, worst `O(m)` |
//' | `"AR1"` | `O(1)` | `O(m)` |
//' | `"AR1Profile"` | `O(1)` | `O(m)` |
//' | `"AR1Focus"` | `O(a_m)`, worst `O(m)` | `O(a_m)`, worst `O(m)` |
//' | `"quantile"` | `O(1)` | `O(1)` |
//' | `"quantileExact"` | expected `O(log m)`, worst `O(m)` | expected `O(log m)`, worst `O(m)` |
//' | `"varCost"` | `O(1)` | `O(1)` |
//' | `"WilcoxonCost"` | `O(m)` | `O(m)` |
//' | `"MedianMoodCost"` | expected `O(log m)`, worst `O(m)` | `O(m)` |
//'
//' For the Gaussian FOCuS test, the standard FOCuS result gives
//' `E[a_m] = O(log m)` under independent continuous noise with a constant or
//' single-change mean, so its complete endpoint update is expected
//' `O(log m)` and worst-case `O(m)`. The recurrence itself is amortized
//' `O(1)`; the logarithmic term comes from maximizing over the active pieces.
//' The Gamma, Gaussian-variance, and AR1Focus implementations use the same
//' active-piece strategy and cost `O(a_m)`, but this package does not claim the
//' Gaussian theorem for those different data models without additional
//' model-specific assumptions.
//'
//' Thus, `"quantile"` and `"varCost"` have genuinely constant-time complete
//' endpoint updates. `"quantileExact"` has an expected logarithmic update.
//' Functional-pruning tests can be fast when `a_m` remains small, but their
//' worst-case endpoint cost is linear.
//'
//' For an AR(1) series generated by `ts_generator(type = "gaussAR1")`, the
//' model is `x[t] = mu[t] + e[t]`, with
//' `e[t] = rho * e[t - 1] + eta[t]` and
//' `eta[t] ~ N(0, sigma2)`. Thus, `rho` controls serial dependence and
//' `sigma2` is the innovation variance, that is, the variance of the new shock
//' after accounting for the previous observation. It is not the marginal
//' variance of the observed series; for a stationary AR(1) noise process, the
//' marginal variance is `sigma2 / (1 - rho^2)`. In `ts_generator()`,
//' `sd_noise` is the innovation standard deviation, so use
//' `sigma2 = sd_noise^2`.
//'
//' `"AR1"` treats `rho` and `sigma2` as fixed and uses the conditional
//' Gaussian likelihood scan given the first observation. `"AR1Profile"` uses
//' the same conditional scan but estimates the innovation variance separately
//' under the no-change and
//' change models. This is useful when the innovation scale is unknown.
//' `"AR1Focus"` is faster because it applies the Gaussian FOCUS calculation
//' to the transformed innovations; it is an approximation and can produce
//' different boundaries near a change point.
//'
//' ## Within-segment cost
//'
//' The validity test fixes the number of segments; `cost` decides between the
//' partitions that attain it, so it controls where the boundaries are placed
//' and not how many there are. Two cost series are available:
//'
//' * `"gaussian"`: the sum of squared deviations of `data[(s + 1):t]` from the
//'   segment mean. This is the right choice for independent observations.
//' * `"ar1"`: the same quantity for the innovations
//'   `z[u] = data[u] - rho * data[u - 1]`. Inside a segment with constant mean
//'   `mu` these have constant mean `(1 - rho) * mu`, so their sum of squared
//'   deviations is the profiled conditional cost of the AR(1) mean model. The
//'   innovation `z[s + 1]` that straddles the boundary is skipped, because its
//'   mean is `mu_new - rho * mu_old` rather than `(1 - rho) * mu_new`; keeping
//'   it would charge every segment for the jump in front of it and pull the
//'   boundary away from the change. Costs are only ever compared across
//'   partitions with the same number of segments `K`, and every such partition
//'   skips exactly `K - 1` innovations, so the comparison remains fair.
//'
//' `"auto"`, the default, uses `"ar1"` for `"AR1"`, `"AR1Profile"`, and
//' `"AR1Focus"`, and `"gaussian"` for every other test. Under serial
//' dependence the Gaussian cost is misspecified, so an AR(1) test combined
//' with the Gaussian cost detects the right number of changes but places them
//' less accurately. Selecting `"ar1"` with a non-AR(1) test is allowed and
//' uses `rho` in the same way.
//'
//' The "subtests" argument selects the validity-based candidate-pruning rules.
//' "right" removes a candidate boundary when its current segment endpoint is
//' invalid; "both" additionally removes that boundary and every older
//' boundary; and "none" retains all candidates. In "right", a valid
//' one-segment candidate is already optimal, so later candidates are deferred
//' until that candidate fails. If the shortcut ends, deferred candidates are
//' replayed through the missed endpoints so that the right-pruning history is
//' preserved; after that, they are tested only at the current endpoint. In
//' "both", candidates are examined from newest to oldest and the scan stops
//' at the first invalid boundary, since the inclusive left rule removes it and
//' all older candidates. In "none", candidates are grouped by their preceding
//' segment count and only the first group containing a valid candidate is
//' needed. Invalid candidates are never used in the current optimum.
//' Pruning is exact only when the selected validity test has the corresponding
//' monotonicity properties. For an arbitrary user-defined rule, use
//' "none" unless those properties have been established.
//'
//' ## Verzelen multiscale validity
//'
//' With `test = "gaussian_mean", subtests = "verzelen"`, a standalone
//' scanner tests every pair of adjacent nonempty intervals contained in the
//' candidate segment. Write \eqn{N} for the full data length and
//' \eqn{C(a,b,c)=\sqrt{(b-a)(c-b)/(c-a)}
//' (\bar y_{a+1:b}-\bar y_{b+1:c})}. The candidate \eqn{(s,t]} is valid iff
//' \deqn{V_s(t)=\max_{s\le a<b<c\le t}\left\{
//' C(a,b,c)^2/(4\sigma^2)-2\log\left(\frac{N(c-a)}{(b-a)(c-b)}\right)
//' \right\}\le\gamma.}
//' Singleton segments have statistic minus infinity and are always valid.
//' The global length \eqn{N} stays fixed when testing shorter segments.
//' Each candidate maintains its current maximum, updated using
//' \eqn{V_s(t)=\max\{V_s(t-1),V_{s+1}(t),H(s,t)\}}, where \eqn{H}
//' scans the splits with outer boundaries \eqn{s,t}. Containment justifies
//' permanent pruning in both directions. FOCUS bounds skip scans that cannot
//' improve the maximum; all splits are checked otherwise. Worst-case total
//' time is cubic in \eqn{N}, with quadratic worst-case state storage.
//'
//' This rule is inspired by the simultaneous CUSUM noise bound on slide 15
//' (PDF page 38) of Nicolas Verzelen's 2023 workshop presentation
//' \emph{Optimal multiple change-point detection and Localization}.
//' It implements the bound with factor 4 and known noise variance `sigma2`.
//' A finite negative `gamma` is allowed for empirical null calibration; such
//' thresholds do not inherit the slide's positive-parameter probability bound.
//' `gamma` is a direct threshold, not a multiplier of `log(N)`.
//'
//' ## MS-l multiscale validity
//'
//' With `test = "gaussian_mean", subtests = "both_msl"`, write \eqn{\ell}
//' for an interval length, \eqn{r} for the shorter side of a split, and
//' \eqn{Z} for its Gaussian two-sample CUSUM. The interval passes when
//' \deqn{\max_{\tau}\left\{|Z_{\tau}|/\sqrt{\sigma^2}
//' -\sqrt{2\log(e\ell/r_{\tau})}\right\}\leq\gamma.}
//' The "both" pruning rule makes a candidate segment valid only when every
//' contained interval passes. Singletons are always valid. `gamma` is the
//' direct threshold \eqn{q}, without a `log(n)` multiplier.
//' The FOCUS maximum bounds the penalized scan; any split that could exceed
//' the threshold is checked exactly.
//' With `subtests = "both_ms_refined"`, the penalty is instead
//' \deqn{\sqrt{2\{1+\log(\ell/r_\tau)+\tfrac12\log(n/\ell)\}},}
//' where \eqn{n} is the full data length. This fixed setting was selected in
//' `simulations/power_gaussian_ms_search` with a separate 1% null calibration
//' and paired F1-curve comparison. For \eqn{n=1000} and \eqn{\sigma^2=1},
//' the calibrated direct threshold is \eqn{q=2.70947832002059}.
//'
//' ## FOCUS iterated-logarithm CUSUM
//'
//' With `subtests = "both_focus_lil"`, each interval of length \eqn{m}
//' uses an envelope of the position penalty of Verzelen et al. (2023,
//' Eq. 12). For a split with side lengths \eqn{l} and \eqn{r=m-l}, put
//' \deqn{h_m(l)=\max\{\min(l,m/l),\min(r,m/r)\},\quad
//' p_m(l)=2\log\log(e h_m(l)).}
//' Its maximum over splits is attained near \eqn{\sqrt m}:
//' \deqn{h_m^*=\max\{\lfloor\sqrt m\rfloor,
//' m/\lceil\sqrt m\rceil\}.}
//' The interval statistic is
//' \deqn{\max_l Z_l^2/\sigma^2-8\log\log(eh_m^*).}
//' FOCUS computes this maximum exactly because the penalty depends only on
//' interval length. This statistic is no larger than the penalized maximum
//' in Proposition 2 with its permitted \eqn{L=2}, so that proposition's
//' single-interval level bound also applies here. A candidate segment is
//' valid when every contained interval has statistic at most `gamma`.
//' The proposition does not directly cover the SVP scan across contained
//' intervals; calibrate `gamma` from the whole-series null maximum for the
//' intended data length and noise model.
//'
//' ## Minimum-arm restricted validity
//'
//' With `subtests = "both_rmin"`, the interval statistic is the same
//' two-sample squared CUSUM \eqn{Z_\tau^2/\sigma^2} used by plain `"both"`,
//' but the scan skips every split whose shorter arm \eqn{r_\tau} is below
//' \eqn{r_{\text{eff}}(m)=\max(r_{\min},\lceil\text{r\_min\_frac}\times
//' m\rceil)} for an interval of current length \eqn{m}:
//' \deqn{\max_{\tau:\,r_\tau\geq r_{\text{eff}}(m)}
//' Z_\tau^2/\sigma^2\leq\gamma.} An interval shorter than
//' \eqn{2\,r_{\text{eff}}(m)} has no eligible split and is always valid.
//' Excluded splits are dropped from the scan rather than down-weighted,
//' unlike `"both_msl"` or `"both_ms_refined"`: this removes the short-scale
//' tests that otherwise dominate the null maximum of a constant threshold,
//' so `gamma` can be calibrated smaller for the same false-detection rate.
//' The absolute floor `r_min` matters most for short intervals; the
//' relative floor `r_min_frac` (zero by default) additionally excludes
//' lopsided splits of long intervals, whose shorter arm is short only
//' relative to the interval. The cost is that a change cannot be resolved
//' within \eqn{r_{\text{eff}}} observations of either end of its segment.
//' `gamma` is a direct threshold on \eqn{Z_\tau^2}, without a `log(n)`
//' multiplier.
//'
//' ## Multiscale CUSUM validity with a simultaneous bound
//'
//' With `subtests = "both_correction"`, a candidate segment \eqn{(s,t]} is valid
//' iff every contained sub-interval passes a length-corrected two-sample
//' CUSUM test:
//' \deqn{T_n(s,t)=\max_{s\le a<u<b\le t}\left\{|Z_u(a,b)|-p_n(u-a,b-u)\right\}
//' \leq\gamma,\qquad p_n(l,r)=\sqrt{2\log\!\left(\frac{n(l+r)}{lr}\right)},}
//' with \eqn{Z_u(a,b)} the standardized two-sample CUSUM comparing
//' \eqn{y_{a+1:u}} to \eqn{y_{u+1:b}} and \eqn{n} the TOTAL series length,
//' fixed for the whole run (never the length of the candidate segment being
//' tested). `gamma` (`q` in the source) is a direct threshold on the
//' \eqn{|Z|} scale, not a `log(n)` multiplier. This is the multiscale "both"
//' construction of `simulations/consistency/SVP_consistency_sharp.tex`,
//' Section 6: the same every-sub-interval scan as `subtests = "verzelen"`
//' (same contrast, same containment argument for exact "both" pruning), but
//' with the penalty \eqn{p_n(l,r)} in place of Verzelen's squared,
//' \eqn{\log(N(c-a)/((b-a)(c-b)))}-penalized statistic. Calibrate `gamma`
//' by simulating the null maximum \eqn{T_n} and taking its empirical
//' \eqn{1-\alpha} quantile (`simulations/power_gaussian/
//' calibrate_consistency_sharp.R`); the paper's simultaneous CUSUM bound
//' keeps this null maximum bounded in probability, so `gamma` need not grow
//' with `n`. Worst-case total time is cubic in the series length, with a
//' FOCUS-based bound skipping scans that cannot improve the running maximum.
//'
//' @param data Numeric vector containing the univariate series. Missing or
//'   non-finite values are not supported.
//' @param gamma Finite scalar validity threshold, positive except when
//'   `subtests` is `"verzelen"`, `"both_msl"`, `"both_ms_refined"`,
//'   `"both_focus_lil"`, `"both_rmin"`, `"both_calibrated"`,
//'   `"both_correction"`, `"both_empirical"`, `"both_lengthonly"`,
//'   `"both_mosum_sqrt"`, or `"both_mosum_adaptive"`. A larger value
//'   generally accepts longer or less homogeneous segments.
//' @param test Character scalar selecting one of the validity tests listed in
//'   Details. Defaults to `"gaussian_mean"`.
//' @param subtests Character scalar selecting the validity-pruning rules:
//'   "both" (default), "right", "none", "verzelen", "both_msl",
//'   "both_ms_refined", "both_focus_lil", "both_rmin", "both_calibrated",
//'   "both_correction", "both_empirical", "both_lengthonly", "both_mosum_sqrt",
//'   "both_mosum", or "both_mosum_adaptive". The last eleven
//'   options select multiscale, minimum-arm-restricted, length-adaptive, or
//'   moving-sum (MOSUM) Gaussian tests; `"both_mosum"` is the MOSUM scan at
//'   its maximal bandwidth (the segment's own midpoint, a single fixed
//'   split), `"both_mosum_sqrt"` the same scan at a smaller,
//'   length-growing bandwidth, and `"both_mosum_adaptive"` the
//'   length-adaptive-threshold version of `"both_mosum"`.
//' @param sigma2 Positive finite innovation variance for fixed-variance AR1
//'   tests (`"AR1"` with `profile_sigma = FALSE`, and `"AR1Focus"`), or
//'   known observation noise variance for `subtests = "verzelen"`,
//'   `"both_msl"`, `"both_ms_refined"`, `"both_focus_lil"`, `"both_rmin"`,
//'   `"both_calibrated"`, `"both_correction"`, `"both_empirical"`,
//'   `"both_lengthonly"`, `"both_mosum_sqrt"`, `"both_mosum"`, and
//'   `"both_mosum_adaptive"` (default 1).
//'   It is the conditional/error variance of the new AR(1) shock, not the
//'   marginal variance of the observed series. It is used when the variance is
//'   fixed. It is ignored by `"AR1Profile"` and by `"AR1"` when
//'   `profile_sigma = TRUE`.
//' @param rho AR(1) coefficient for the three AR1 tests. It must be finite and
//'   strictly between -1 and 1. If it is `NA`, the value is estimated robustly
//'   from the full series using [AR1_rho()].
//' @param profile_sigma Logical; if `TRUE`, estimate the innovation variance
//'   separately under the no-change and one-change AR(1) models when
//'   `test = "AR1"`. This removes the need to know the innovation scale and
//'   uses the resulting residual sums of squares in the likelihood-ratio
//'   statistic. Using `test = "AR1Profile"` has the same effect. The argument
//'   is ignored by other tests.
//' @param quantile Quantile level used by `"quantile"` and
//'   `"quantileExact"`. It is ignored by other tests.
//' @param cost Character scalar selecting the within-segment cost used to rank
//'   partitions with the same number of segments: `"auto"` (default),
//'   `"gaussian"`, or `"ar1"`. See Details. `"ar1"` uses `rho`, estimating it
//'   from the full series with [AR1_rho()] when it is `NA`.
//' @param r_min Positive integer absolute minimum arm length for
//'   `subtests = "both_rmin"`. Ignored by every other `subtests` value.
//' @param r_min_frac Finite scalar in `[0, 0.5)`, the minimum arm length for
//'   `subtests = "both_rmin"` as a fraction of the interval's current
//'   length, combined with `r_min` by taking the larger of the two floors.
//'   Zero (default) disables the relative floor. Ignored by every other
//'   `subtests` value.
//' @param abc_A,abc_B,abc_C Coefficients of the length-dependent penalty
//'   `abc_A * log(m) + abc_B * log(log(m)) + abc_C` subtracted from the
//'   Gaussian statistic of a candidate segment of length `m` when
//'   `subtests = "both_calibrated"`. `NA` (default) uses the built-in
//'   calibrated values 1.682969, 6.208626 and 4.939248. Ignored by every
//'   other `subtests` value.
//'
//' @return A list with "changepoints" (the inclusive end of every segment,
//'   including "length(data)"), "lastIndexSet" (zero-based candidate boundaries
//'   remaining at termination, in decreasing order), "nb" (the number of
//'   active candidates at the beginning of each endpoint iteration, before
//'   pruning), "costQ" (always "NULL"), and "R". Row "t" of matrix "R"
//'   stores the best cumulative cost on the scale selected by `cost`, the
//'   number of segments, and the previous boundary for `data[1:t]`. The three
//'   AR(1) tests also return "rho" and "sigma2".
//'
//' @examples
//' # Gaussian mean: FOCuS test for independent Gaussian observations.
//' set.seed(1)
//' gaussian_data <- ts_generator(
//'   chpts = c(20, 40, 60), parameters = c(0, 2, -1),
//'   sd_noise = 1, type = "gauss"
//' )
//' SVP(gaussian_data, gamma = 4 * log(length(gaussian_data)),
//'     test = "gaussian_mean")$changepoints
//'
//' # Gamma rate: positive exponential observations with changing rates.
//' gamma_data <- ts_generator(
//'   chpts = c(20, 40, 60), parameters = c(1, 4, 2),
//'   type = "exp"
//' )
//' SVP(gamma_data, gamma = 2 * log(length(gamma_data)),
//'     test = "gamma_rate")$changepoints
//'
//' # Gaussian variance: parameters are segment standard deviations.
//' variance_data <- ts_generator(
//'   chpts = c(20, 40, 60), parameters = c(0.5, 1.5, 0.75),
//'   type = "variance"
//' )
//' SVP(variance_data, gamma = 2 * log(length(variance_data)),
//'     test = "gaussian_variance")$changepoints
//'
//' # Quantile and exact quantile tests: robust tests for changes in spread.
//' quantile_data <- ts_generator(
//'   chpts = c(20, 40, 60), parameters = c(0.5, 1.5, 0.75),
//'   type = "variance"
//' )
//' SVP(quantile_data, gamma = 2, quantile = 0.1,
//'     test = "quantile")$changepoints
//' SVP(quantile_data, gamma = 2, quantile = 0.1,
//'     test = "quantileExact")$changepoints
//'
//' # Robust variance, Wilcoxon, and Median-Mood tests.
//' robust_data <- ts_generator(
//'   chpts = c(20, 40, 60), parameters = c(0, 2, -1),
//'   sd_noise = 1, type = "gauss"
//' )
//' SVP(robust_data, gamma = 2, test = "varCost")$changepoints
//' SVP(robust_data, gamma = 2 * log(length(robust_data)),
//'     test = "WilcoxonCost")$changepoints
//' SVP(robust_data, gamma = 2 * log(length(robust_data)),
//'     test = "MedianMoodCost")$changepoints
//'
//' # Exact AR(1), with the known innovation variance.
//' ar1_data <- ts_generator(
//'   chpts = c(20, 40, 60), parameters = c(0, 2, -1),
//'   sd_noise = 0.8, rho = 0.7, type = "gaussAR1"
//' )
//' SVP(ar1_data, gamma = 4 * log(length(ar1_data)), test = "AR1",
//'     rho = 0.7, sigma2 = 0.8^2)$changepoints
//'
//' # Exact AR(1) with the innovation variance profiled out.
//' SVP(ar1_data, gamma = 2 * log(length(ar1_data)),
//'     test = "AR1Profile", rho = 0.7, sigma2 = 1)$changepoints
//'
//' # Faster approximate AR(1) FOCUS test on the transformed innovations.
//' SVP(ar1_data, gamma = 2 * log(length(ar1_data)),
//'     test = "AR1Focus", rho = 0.7, sigma2 = 0.8^2)$changepoints
//'
//' # The same test ranking partitions with the Gaussian cost instead, which
//' # keeps the number of segments but can move the boundaries.
//' SVP(ar1_data, gamma = 2 * log(length(ar1_data)),
//'     test = "AR1Focus", rho = 0.7, sigma2 = 0.8^2,
//'     cost = "gaussian")$changepoints
//'
//' @seealso [svp0()] for arbitrary R validity functions, [AR1_rho()], and
//'   [AR1_single_change()].
//' @export
// [[Rcpp::export]]
List SVP(std::vector<double> data,
         double gamma,
         std::string test = "gaussian_mean",
         std::string subtests = "both",
         double sigma2 = 1.0,
         double rho = NA_REAL,
         bool profile_sigma = false,
         double quantile = 0.01,
         std::string cost = "auto",
         int r_min = 1,
         double r_min_frac = 0.0,
         double abc_A = NA_REAL,
         double abc_B = NA_REAL,
         double abc_C = NA_REAL)
{
  if (data.empty()) {
    stop("'data' must contain at least one observation.");
  }

  if (!std::all_of(data.begin(), data.end(),
                   [](double value) { return std::isfinite(value); })) {
    stop("'data' must contain only finite observations.");
  }

  const bool multiscale = subtests == "verzelen" ||
    subtests == "both_msl" || subtests == "both_ms_refined" ||
    subtests == "both_focus_lil" || subtests == "both_rmin" ||
    subtests == "both_calibrated" ||
    subtests == "both_correction" || subtests == "both_empirical" ||
    subtests == "both_lengthonly" || subtests == "both_mosum_sqrt" ||
    subtests == "both_mosum_adaptive";
  if (!std::isfinite(gamma) || (!multiscale && gamma <= 0.0)) {
    if (multiscale)
      stop("gamma must be finite");
    stop("gamma must be finite and positive");
  }

  const TestKind test_kind = parse_test_kind(test);
  const CostKind cost_kind = parse_cost_kind(cost);

  if (subtests != "both" && subtests != "right" && subtests != "none" &&
      subtests != "verzelen" && subtests != "both_msl" &&
      subtests != "both_ms_refined" && subtests != "both_focus_lil" &&
      subtests != "both_rmin" &&
      subtests != "both_calibrated" && subtests != "both_correction" &&
      subtests != "both_empirical" && subtests != "both_lengthonly" &&
      subtests != "both_mosum_sqrt" && subtests != "both_mosum" &&
      subtests != "both_mosum_adaptive") {
    stop("'subtests' must be one of 'both', 'right', 'none', 'verzelen', 'both_msl', 'both_ms_refined', 'both_focus_lil', 'both_rmin', 'both_calibrated', 'both_correction', 'both_empirical', 'both_lengthonly', 'both_mosum_sqrt', 'both_mosum', or 'both_mosum_adaptive'.");
  }
  if (multiscale) {
    if (test_kind != TestKind::GaussianMean)
      stop("multiscale subtests require test = 'gaussian_mean'");
    if (!std::isfinite(sigma2) || sigma2 <= 0.0)
      stop("sigma2 must be finite and positive");
  }
  if (subtests == "both_rmin" && r_min < 1) {
    stop("r_min must be at least one");
  }
  if (subtests == "both_rmin" &&
      (!std::isfinite(r_min_frac) || r_min_frac < 0.0 || r_min_frac >= 0.5)) {
    stop("r_min_frac must be finite and in [0, 0.5)");
  }

  if ((test_kind == TestKind::Quantile ||
       test_kind == TestKind::QuantileExact) &&
      (!std::isfinite(quantile) || quantile <= 0.0 || quantile > 0.5)) {
    stop("quantile must be finite and in (0, 0.5]");
  }

  if (test_kind == TestKind::GammaRate &&
      !std::all_of(data.begin(), data.end(),
                   [](double value) { return value > 0.0; })) {
    stop("'data' must be strictly positive for test 'gamma_rate'.");
  }

  const bool ar1_test = test_kind == TestKind::AR1 ||
    test_kind == TestKind::AR1Profile ||
    test_kind == TestKind::AR1Focus;
  const bool ar1_cost = cost_kind == CostKind::AR1 ||
    (cost_kind == CostKind::Auto && ar1_test);

  // The AR(1) cost and the three AR(1) tests share one rho, so it is resolved
  // once here. It is estimated from the full series when it is not supplied.
  double rho_used = rho;
  if (ar1_test || ar1_cost) {
    if (NumericVector::is_na(rho_used)) {
      rho_used = svp_detail::robust_ar1_rho(data);
    }
    if (!std::isfinite(rho_used) || std::fabs(rho_used) >= 1.0) {
      stop("rho must be finite and strictly between -1 and 1");
    }
  }

  const svp_detail::SegmentCost segment_cost = ar1_cost ?
    svp_detail::SegmentCost::ar1(data, rho_used) :
    svp_detail::SegmentCost::gaussian(data);

  switch (test_kind) {
  case TestKind::GaussianMean:
    if (subtests == "verzelen")
      return svp_detail::svp_verzelen_impl(data, gamma, sigma2, segment_cost);
    if (subtests == "both_correction")
      return svp_detail::svp_consistency_sharp_impl(data, gamma, sigma2, segment_cost);
    if (subtests == "both_empirical")
      return svp_detail::svp_empirical_both_impl(data, gamma, sigma2, segment_cost);
    if (subtests == "both_msl")
      return svp_detail::svp_impl<GaussianMeanMSL>(
        data, gamma, "both", segment_cost, sigma2
      );
    if (subtests == "both_ms_refined")
      return svp_detail::svp_impl<GaussianMeanMSFamily>(
        data, gamma, "both", segment_cost, sigma2, data.size(), 1.0, 0.5
      );
    if (subtests == "both_focus_lil")
      return svp_detail::svp_impl<GaussianMeanFocusLIL>(
        data, gamma, "both", segment_cost, sigma2
      );
    if (subtests == "both_rmin")
      return svp_detail::svp_impl<GaussianMeanRMin>(
        data, gamma, "both", segment_cost, sigma2,
        static_cast<std::size_t>(r_min), r_min_frac
      );
    if (subtests == "both_calibrated")
      return svp_detail::svp_impl<GaussianMeanBothCalibrated>(
        data, gamma, "both", segment_cost, sigma2,
        NumericVector::is_na(abc_A) ? GaussianMeanBothCalibrated::kDefaultA : abc_A,
        NumericVector::is_na(abc_B) ? GaussianMeanBothCalibrated::kDefaultB : abc_B,
        NumericVector::is_na(abc_C) ? GaussianMeanBothCalibrated::kDefaultC : abc_C
      );
    if (subtests == "both_lengthonly")
      return svp_detail::svp_impl<GaussianMeanLengthOnly>(
        data, gamma, "both", segment_cost, sigma2,
        static_cast<double>(data.size())
      );
    if (subtests == "both_mosum_sqrt")
      return svp_detail::svp_impl<GaussianMeanMosumSqrt>(
        data, gamma, "both", segment_cost, sigma2
      );
    if (subtests == "both_mosum")
      return svp_detail::svp_impl<GaussianMeanMosum>(
        data, gamma, "both", segment_cost, sigma2
      );
    if (subtests == "both_mosum_adaptive")
      return svp_detail::svp_impl<GaussianMeanMosumAdaptive>(
        data, gamma, "both", segment_cost, sigma2
      );
    return svp_detail::svp_impl<GaussianMean>(data, gamma, subtests,
                                              segment_cost);
  case TestKind::GammaRate:
    return svp_detail::svp_impl<GammaRate>(data, gamma, subtests,
                                          segment_cost);
  case TestKind::GaussianVariance:
    return svp_detail::svp_impl<GaussianVariance>(data, gamma, subtests,
                                                 segment_cost);
  case TestKind::QuantileExact:
    return svp_detail::svp_impl<QuantileCostExact>(
      data, gamma, subtests, segment_cost, quantile
    );
  case TestKind::Quantile:
    return svp_detail::svp_impl<QuantileCost>(
      data, gamma, subtests, segment_cost, quantile
    );
  case TestKind::VarCost:
    return svp_detail::svp_impl<VarianceCost>(data, gamma, subtests,
                                             segment_cost);
  case TestKind::WilcoxonCost:
    return svp_detail::svp_impl<WilcoxonCost>(data, gamma, subtests,
                                             segment_cost);
  case TestKind::MedianMoodCost:
    return svp_detail::svp_impl<MedianMoodCost>(data, gamma, subtests,
                                               segment_cost);
  case TestKind::AR1:
  case TestKind::AR1Profile:
  case TestKind::AR1Focus: {
    const bool profiles_sigma = test_kind == TestKind::AR1Profile ||
      (test_kind == TestKind::AR1 && profile_sigma);
    if (!profiles_sigma &&
        (!std::isfinite(sigma2) || sigma2 <= 0.0)) {
      stop("sigma2 must be finite and positive");
    }
    if (test_kind == TestKind::AR1Focus) {
      List result = svp_detail::svp_impl<AR1FocusMeanChange>(
        data,
        gamma,
        subtests,
        segment_cost,
        rho_used,
        sigma2
      );
      result["rho"] = rho_used;
      result["sigma2"] = sigma2;
      return result;
    }
    List result = svp_detail::svp_impl<AR1ExactMeanChange>(
      data,
      gamma,
      subtests,
      segment_cost,
      rho_used,
      sigma2,
      profiles_sigma
    );
    result["rho"] = rho_used;
    result["sigma2"] = sigma2;
    return result;
  }
  }

  stop("Internal error while selecting the SVP test.");
}
