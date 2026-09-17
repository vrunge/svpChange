% Generated from the documentation in src/svp.cpp
\name{SVP}
\alias{SVP}
\title{Smallest Valid Partitioning with Incremental Validity Tests}
\usage{
SVP(
  data,
  gamma,
  test = "gaussian_mean",
  subtests = "both",
  sigma2 = 1,
  rho = NA_real_,
  profile_sigma = FALSE,
  quantile = 0.01,
  cost = "auto"
)
}
\arguments{
\item{data}{Numeric vector containing the univariate series. Missing and
non-finite values are not supported.}

\item{gamma}{Positive finite scalar validity threshold. Larger values
generally accept longer or less homogeneous segments.}

\item{test}{Character scalar selecting a validity test; see Details. Defaults
to \code{"gaussian_mean"}.}

\item{subtests}{Character scalar selecting the validity-pruning rules:
\code{"both"} (default), \code{"right"}, or \code{"none"}.}

\item{sigma2}{Positive finite innovation variance for fixed-variance AR1
tests (\code{"AR1"} with \code{profile_sigma = FALSE}, and
\code{"AR1Focus"}). It is the conditional/error variance of the new AR(1)
shock, not the marginal variance of the observed series. It is ignored by
\code{"AR1Profile"} and by \code{"AR1"} when
\code{profile_sigma = TRUE}.}

\item{rho}{AR(1) coefficient for the AR1 tests. It must be finite and strictly
  between -1 and 1. If it is \code{NA}, it is estimated robustly from the full
  series using \code{AR1_rho()}.}

\item{profile_sigma}{If \code{TRUE}, estimate the innovation variance
separately under the no-change and one-change AR(1) models when
\code{test = "AR1"}. This removes the need to know the innovation scale and
uses the resulting residual sums of squares in the likelihood-ratio
statistic. Using \code{test = "AR1Profile"} has the same effect. The argument
is ignored by other tests.}

\item{quantile}{Quantile level used by \code{"quantile"} and
\code{"quantileExact"}; ignored by other tests.}

\item{cost}{Character scalar selecting the within-segment cost used to rank
partitions with the same number of segments: \code{"auto"} (default),
\code{"gaussian"}, or \code{"ar1"}. See Details. \code{"ar1"} uses
\code{rho}, estimating it from the full series with \code{AR1_rho()} when it
is \code{NA}.}
}
\value{
A list containing:
\describe{
  \item{changepoints}{Inclusive end of every segment, including
  \code{length(data)}.}
  \item{lastIndexSet}{Candidate boundaries remaining at termination.}
  \item{nb}{Number of active candidate boundaries at the beginning of each
  endpoint iteration, before pruning.}
  \item{costQ}{Currently \code{NULL}.}
  \item{R}{Dynamic-programming matrix. Row \code{t} contains the best
  cumulative cost on the scale selected by \code{cost}, the number of
  segments, and the previous boundary for \code{data[1:t]}.}
  \item{rho, sigma2}{Returned by the three AR(1) tests only.}
}
}
\description{
Segments a univariate series into the smallest number of segments that pass a
selected validity test. A candidate segment is evaluated only at its current
endpoint; its incremental state receives intervening observations, but earlier
endpoints are not tested again. With \code{subtests = "none"}, an invalid
endpoint does not remove the candidate and it can recover later. With
\code{subtests = "right"} or \code{subtests = "both"}, it is removed
immediately. With \code{subtests = "both"}, an invalid boundary also removes
that boundary and every older boundary under the inclusive left-pruning rule.
Among partitions with the same number of
segments, it chooses the one with the smallest within-segment cost, selected
with \code{cost}. Validity statistics are maintained incrementally in C++.
}
\details{
The available validity tests are:
\itemize{
  \item \code{"gaussian_mean"}: Gaussian FOCUS likelihood-ratio statistic
  evaluated at the current segment endpoint for a change in mean. Use for
  independent Gaussian data with constant variance. To require every prefix
  to be valid, use \code{svp0()} with \code{valid_FOCUS()}.
  \item \code{"gamma_rate"}: likelihood-ratio test for a rate change in
  positive Gamma data (unit shape in the current implementation).
  \item \code{"gaussian_variance"}: the Gamma-rate test applied to squared
  observations, for variance changes around a known zero mean.
  \item \code{"AR1"}: conditional fixed-\code{rho} Gaussian likelihood-ratio
  scan for a change in the marginal mean of an AR(1) series, conditional on the
  first observation. In the AR(1) model, an
  innovation is the new random shock after accounting for the previous
  observation and the AR(1) mean structure; \code{sigma2} is the variance of
  this shock, not the marginal variance of the observations.
  \item \code{"AR1Profile"}: the same conditional AR(1) scan, but profiling out the
  innovation variance. The variance is estimated separately under the
  no-change and change models from their residual sums of squares, which is
  useful when the innovation scale is unknown.
  \item \code{"AR1Focus"}: faster approximation applying Gaussian FOCUS to
  innovations \code{x[t] - rho * x[t - 1]}. Prefer \code{"AR1"} when exact
  boundary treatment matters.
  \item \code{"quantile"} and \code{"quantileExact"}: quantile-based robust
  tests, using \code{quantile} as their quantile level.
  \item \code{"varCost"}: robust variance cost test.
  \item \code{"WilcoxonCost"}: Wilcoxon rank-based cost test.
  \item \code{"MedianMoodCost"}: Median-Mood rank-based cost test.
}

\subsection{Per-candidate update complexity}{
Let \eqn{m} be the current length of one candidate segment and let
\eqn{a_m} be the number of active functional-pruning pieces, with
\eqn{a_m \le m}. A state update incorporates one new observation. A complete
endpoint check also computes the statistic needed to decide validity. These
are per-candidate costs; total SVP time also depends on the number of
candidates examined under the selected \code{subtests} mode.

\tabular{lll}{
Test \tab State update \tab Update and endpoint check \cr
\code{"gaussian_mean"} \tab amortized \eqn{O(1)}, worst \eqn{O(m)}
  \tab expected \eqn{O(\log m)}, worst \eqn{O(m)} \cr
\code{"gamma_rate"} \tab \eqn{O(a_m)}, worst \eqn{O(m)}
  \tab \eqn{O(a_m)}, worst \eqn{O(m)} \cr
\code{"gaussian_variance"} \tab \eqn{O(a_m)}, worst \eqn{O(m)}
  \tab \eqn{O(a_m)}, worst \eqn{O(m)} \cr
\code{"AR1"} \tab \eqn{O(1)} \tab \eqn{O(m)} \cr
\code{"AR1Profile"} \tab \eqn{O(1)} \tab \eqn{O(m)} \cr
\code{"AR1Focus"} \tab \eqn{O(a_m)}, worst \eqn{O(m)}
  \tab \eqn{O(a_m)}, worst \eqn{O(m)} \cr
\code{"quantile"} \tab \eqn{O(1)} \tab \eqn{O(1)} \cr
\code{"quantileExact"} \tab expected \eqn{O(\log m)}, worst \eqn{O(m)}
  \tab expected \eqn{O(\log m)}, worst \eqn{O(m)} \cr
\code{"varCost"} \tab \eqn{O(1)} \tab \eqn{O(1)} \cr
\code{"WilcoxonCost"} \tab \eqn{O(m)} \tab \eqn{O(m)} \cr
\code{"MedianMoodCost"} \tab expected \eqn{O(\log m)}, worst \eqn{O(m)}
  \tab \eqn{O(m)} \cr
}

For the Gaussian FOCuS test, the standard FOCuS result gives
\eqn{E[a_m] = O(\log m)} under independent continuous noise with a constant
or single-change mean. Its complete endpoint update is therefore expected
\eqn{O(\log m)} and worst-case \eqn{O(m)}. The recurrence itself is amortized
\eqn{O(1)}; the logarithmic term comes from maximizing over the active pieces.
The Gamma, Gaussian-variance, and AR1Focus implementations use the same
active-piece strategy and cost \eqn{O(a_m)}, but this package does not claim
the Gaussian theorem for those different data models without additional
model-specific assumptions.

The genuinely constant-time complete endpoint updates are \code{"quantile"}
and \code{"varCost"}. \code{"quantileExact"} has expected logarithmic cost.
Functional-pruning tests can be fast when \eqn{a_m} remains small, but their
worst-case endpoint cost is linear.
}

For an AR(1) series generated by \code{ts_generator(type = "gaussAR1")}, the
model is \code{x[t] = mu[t] + e[t]}, with
\code{e[t] = rho * e[t - 1] + eta[t]} and
\code{eta[t] ~ N(0, sigma2)}. Thus, \code{rho} controls serial dependence and
\code{sigma2} is the innovation variance, that is, the variance of the new shock
after accounting for the previous observation. It is not the marginal
variance of the observed series; for a stationary AR(1) noise process, the
marginal variance is \code{sigma2 / (1 - rho^2)}. In \code{ts_generator()},
\code{sd_noise} is the innovation standard deviation, so use
\code{sigma2 = sd_noise^2}.

\code{"AR1"} treats \code{rho} and \code{sigma2} as fixed and uses the exact
conditional Gaussian likelihood scan given the first observation.
\code{"AR1Profile"} uses the same conditional scan but estimates the innovation variance separately under the no-change and
change models. This is useful when the innovation scale is unknown.
\code{"AR1Focus"} is faster because it applies the Gaussian FOCUS calculation
to the transformed innovations; it is an approximation and can produce
different boundaries near a change point.

\subsection{Within-segment cost}{
The validity test fixes the number of segments; \code{cost} decides between
the partitions that attain it, so it controls where the boundaries are placed
and not how many there are. Two cost series are available:
\itemize{
  \item \code{"gaussian"}: the sum of squared deviations of
  \code{data[(s + 1):t]} from the segment mean. This is the right choice for
  independent observations.
  \item \code{"ar1"}: the same quantity for the innovations
  \code{z[u] = data[u] - rho * data[u - 1]}. Inside a segment with constant
  mean \code{mu} these have constant mean \code{(1 - rho) * mu}, so their sum
  of squared deviations is the profiled conditional cost of the AR(1) mean
  model. The innovation \code{z[s + 1]} that straddles the boundary is
  skipped, because its mean is \code{mu_new - rho * mu_old} rather than
  \code{(1 - rho) * mu_new}; keeping it would charge every segment for the
  jump in front of it and pull the boundary away from the change. Costs are
  only ever compared across partitions with the same number of segments
  \code{K}, and every such partition skips exactly \code{K - 1} innovations,
  so the comparison remains fair.
}

\code{"auto"}, the default, uses \code{"ar1"} for \code{"AR1"},
\code{"AR1Profile"}, and \code{"AR1Focus"}, and \code{"gaussian"} for every
other test. Under serial dependence the Gaussian cost is misspecified, so an
AR(1) test combined with the Gaussian cost detects the right number of changes
but places them less accurately. Selecting \code{"ar1"} with a non-AR(1) test
is allowed and uses \code{rho} in the same way.
}

The \code{subtests} argument selects validity-based candidate pruning.
\code{"right"} removes a candidate when its current endpoint is invalid;
when the one-segment candidate is valid, later candidates are deferred until
it fails. If the shortcut ends, deferred candidates are replayed through the
missed endpoints so that the right-pruning history is preserved; after that,
they are tested only at the current endpoint. \code{"both"} additionally
removes that boundary and every older boundary; candidates are examined from
newest to oldest and the scan stops at the first invalid boundary.
\code{"none"} retains all candidates and groups them by their predecessor
segment count, stopping after the first group containing a valid candidate.
Invalid candidates are never used in the current optimum.
Pruning is exact only when the selected validity statistic has the
corresponding monotonicity properties.

The available nonparametric and robust cost tests are \code{"quantile"},
\code{"quantileExact"}, \code{"varCost"}, \code{"WilcoxonCost"}, and
\code{"MedianMoodCost"}. They are selected through the same \code{test}
argument and use the same pruning options as the Gaussian tests. The
\code{nb} output records the number of active candidates at the beginning of
each endpoint iteration, before pruning.
}
\examples{
# Gaussian mean: FOCuS test for independent Gaussian observations.
set.seed(1)
gaussian_data <- ts_generator(
  chpts = c(20, 40, 60), parameters = c(0, 2, -1),
  sd_noise = 1, type = "gauss"
)
SVP(gaussian_data, gamma = 2 * log(length(gaussian_data)),
    test = "gaussian_mean")$changepoints

# Gamma rate: positive exponential observations with changing rates.
gamma_data <- ts_generator(
  chpts = c(20, 40, 60), parameters = c(1, 4, 2),
  type = "exp"
)
SVP(gamma_data, gamma = 2 * log(length(gamma_data)),
    test = "gamma_rate")$changepoints

# Gaussian variance: parameters are segment standard deviations.
variance_data <- ts_generator(
  chpts = c(20, 40, 60), parameters = c(0.5, 1.5, 0.75),
  type = "variance"
)
SVP(variance_data, gamma = 2 * log(length(variance_data)),
    test = "gaussian_variance")$changepoints

# Quantile and exact quantile tests: robust tests for changes in spread.
quantile_data <- ts_generator(
  chpts = c(20, 40, 60), parameters = c(0.5, 1.5, 0.75),
  type = "variance"
)
SVP(quantile_data, gamma = 2, quantile = 0.1,
    test = "quantile")$changepoints
SVP(quantile_data, gamma = 2, quantile = 0.1,
    test = "quantileExact")$changepoints

# Robust variance, Wilcoxon, and Median-Mood tests.
robust_data <- ts_generator(
  chpts = c(20, 40, 60), parameters = c(0, 2, -1),
  sd_noise = 1, type = "gauss"
)
SVP(robust_data, gamma = 2, test = "varCost")$changepoints
SVP(robust_data, gamma = 2 * log(length(robust_data)),
    test = "WilcoxonCost")$changepoints
SVP(robust_data, gamma = 2 * log(length(robust_data)),
    test = "MedianMoodCost")$changepoints

# Exact AR(1), with the known innovation variance.
ar1_data <- ts_generator(
  chpts = c(20, 40, 60), parameters = c(0, 2, -1),
  sd_noise = 0.8, rho = 0.7, type = "gaussAR1"
)
SVP(ar1_data, gamma = 2 * log(length(ar1_data)), test = "AR1",
    rho = 0.7, sigma2 = 0.8^2)$changepoints

# Exact AR(1) with the innovation variance profiled out.
SVP(ar1_data, gamma = 2 * log(length(ar1_data)),
    test = "AR1Profile", rho = 0.7, sigma2 = 1)$changepoints

# Faster approximate AR(1) FOCUS test on the transformed innovations.
SVP(ar1_data, gamma = 2 * log(length(ar1_data)),
    test = "AR1Focus", rho = 0.7, sigma2 = 0.8^2)$changepoints

# The same test ranking partitions with the Gaussian cost instead, which keeps
# the number of segments but can move the boundaries.
SVP(ar1_data, gamma = 2 * log(length(ar1_data)),
    test = "AR1Focus", rho = 0.7, sigma2 = 0.8^2,
    cost = "gaussian")$changepoints
}
\seealso{
\code{\link{svp0}}, \code{\link{AR1_rho}}, \code{\link{AR1_single_change}}
}
