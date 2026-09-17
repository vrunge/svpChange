# svpChange

`svpChange` implements Smallest Valid Partitioning (SVP) for multiple
change-point detection in univariate time series. SVP combines dynamic
programming with a local validity constraint: a segment is retained only when
its selected validity statistic stays below a user-defined threshold.

The computational validity tests are implemented in C++ through Rcpp.

## Installation

```r
install.packages("remotes")
remotes::install_github("vrunge/svpChange")
```

For local development:

```r
remotes::install_local("path/to/svpChange")
```

## Basic use

```r
library(svpChange)
set.seed(1)
n <- 120
y <- ts_generator(
  chpts = c(40, 80, 120), parameters = c(0, 2, -1),
  sd_noise = 1, type = "gauss"
)

fit <- SVP(y, gamma = 1.5 * log(length(y)), test = "gaussian_mean")
fit$changepoints
```

`changepoints` contains the inclusive end index of every segment, including
the final observation. `R` contains the dynamic-programming cost, segment
count, and previous boundary for each time point.

## Built-in validity tests

`SVP()` accepts these values of `test`:

| Test | Use |
| --- | --- |
| `gaussian_mean` | Gaussian FOCUS mean-change test |
| `gamma_rate` | Change in rate for positive Gamma observations |
| `gaussian_variance` | Variance changes through squared observations |
| `quantile`, `quantileExact` | Quantile-based robust tests |
| `varCost` | Robust variance cost test |
| `WilcoxonCost` | Wilcoxon rank-based test |
| `MedianMoodCost` | Median-Mood rank-based test |
| `AR1` | Exact fixed-`rho` AR(1) mean-change test |
| `AR1Profile` | AR(1) test with profiled innovation variance |
| `AR1Focus` | Faster innovation-based AR(1) approximation |

For quantile tests:

```r
SVP(y, gamma = 10, test = "quantile", quantile = 0.05)
```

For AR(1) data:

```r
fit_ar1 <- SVP(y, gamma = 10, test = "AR1", rho = 0.7, sigma2 = 1)
```

If `rho` is omitted, `SVP()` estimates it robustly using `AR1_rho()`.

## Pruning and multiscale thresholds

The `subtests` argument controls the candidate set:

```r
SVP(y, gamma = 10, test = "gaussian_mean",
    subtests = "both")
```

The values are `"none"`, `"right"`, and `"both"`. Use `"none"`
for arbitrary validity rules unless the pruning assumptions have been
established for the selected test.

## User-defined validity tests

Use `svp0()` when the validity rule is an R function:

```r
my_test <- function(segment, gamma) max(segment) - min(segment) <= gamma
fit <- svp0(y, gamma = 5, test = my_test)
```

The package also provides `valid_FOCUS()`, `valid_AR1()`, `valid_SSE()`,
`valid_RANGE()`, `valid_RANGE_SLACK()`, `valid_QUANTILE()`, `valid_SCALE()`,
and `valid_OP()` for `svp0()`.

## Other algorithms and utilities

`OP()`, `PELT()`, and `SN()` are available for algorithm comparisons.
`AR1_rho()` and `AR1_single_change()` provide AR(1) diagnostics, and
`ts_generator()` generates simulation signals.

## Simulations

The `simulations/` directory contains paper-style experiments:

```text
power_gaussian/  Gaussian power studies
power_ar1/       AR(1) power studies
power_robust/    Heavy-tailed power studies
time_gaussian/   Gaussian runtime studies
time_ar1/        AR(1) runtime studies
time_robust/     Robust runtime studies
other_simus/     Supporting, historical, and application studies
```

Each primary folder contains a `README.md`. Large Monte Carlo studies are opt-in; see
`simulations/README.md` for the run convention.

## Testing

From the package root:

```r
devtools::test()
```

The tests cover Gaussian, robust, quantile, AR(1), pruning, and
algorithm-equivalence behavior.

## License

GPL-3. See `DESCRIPTION` for authorship and package metadata.
