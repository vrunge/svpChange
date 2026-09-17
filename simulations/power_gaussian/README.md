# Gaussian power study

The maintained Section 5.1 study uses `n = 1000`, the four paper signals,
jump sizes from 0.1 to 2, 100 replications, independent standard Gaussian
noise, and a matching tolerance of `round(0.0025*n)`. The deterministic signal
patterns and Gaussian noise generator match the preserved paper simulation.

The four plotted methods are PELT with penalty `2*log(n)`, SVP Gaussian /
right / c = 2.000, and the calibrated `right` and `both` modes. The calibrated
constants are read from `gaussian_mode_calibration.csv`. They match the PELT
null false-positive probability on 20,000 independent calibration series and
are checked on 10,000 separate validation series. The `none` mode is calibrated
by the same protocol but is reserved for the timing comparison.

Changepoint scores exclude the mandatory terminal boundary `n`. Matching is
one-to-one within the stated tolerance. Consequently, an unchanged series has
F1 equal to one only when no internal changepoint is reported.

Run from the package root:

```r
Sys.setenv(SVP_GAUSSIAN_CALIBRATION_WORKERS = "8")
source("simulations/power_gaussian/calibrate_modes.R")
run_gaussian_mode_calibration(workers = 8L)

Sys.setenv(SVP_RUN_SIMULATIONS = "false")
source("simulations/power_gaussian/run_power.R")
run_and_save_gaussian(workers = 8L)
source("simulations/power_gaussian/export_paper_figures.R")
export_gaussian_paper_figures()
```

`compare_gaussian_refresh.R` separates the scoring correction from the
constant recalibration. Set `SVP_COMPARE_GAUSSIAN=true` when sourcing it to
write `gaussian_performance_decomposition.csv` and
`gaussian_performance_changes.csv`.

The historical current-package reproduction is retained in
`reproduced_v020_matched/`. It is not the source of the refreshed figures.
