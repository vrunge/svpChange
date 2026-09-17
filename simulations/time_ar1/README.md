# AR(1) runtime study

`run_time.R` benchmarks the AR(1) power-study competitors under stationary
Gaussian AR(1) noise with rho = 0.8:

- `PELT AR1 approximate` with its calibrated innovation-PELT penalty;
- `PELT inflated` with the fixed heuristic inflation formula;
- calibrated `DeCAFS AR1`;
- calibrated `SVP AR1Focus + AR1 cost / right`.

It varies sequence length under no change and the number of changes at fixed
length. Calibration constants are read from `power_ar1/ar1_calibration.rds`.

The sequence-length grid contains 24 values from 256 through 32,768; the
fixed-length experiment uses `n = 2000`, jump size 10, and `k = 0:99`.
There are 64 replicates per method and configuration. Fast fits are
adaptively batched for approximately 0.05 seconds before dividing by the
batch size.
The preferred paper output is `ar1_time_time_pair.pdf`: a side-by-side
two-panel figure with one shared legend above the panels. The individual
outputs remain `ar1_time_time_vs_n.pdf` and
`ar1_time_time_vs_detected.pdf` for compatibility. The first plot uses the median elapsed time
under the null; the second uses the median at each exact detected segment count
from the fixed-`n` experiment, with bootstrap 95% intervals, displaying
detected segments up to 100. The old
`time_vs_k` file, if present, is historical and is not regenerated.

Because the two PELT fits are much faster than the SVP and DeCAFS fits, their
detected-segment timings use a separate rerun with 16 replicates and a
0.2-second adaptive timing batch. The DeCAFS and SVP timings retain the main
64-replicate, 0.05-second-batch run.

Run from the package root:

```r
Sys.setenv(SVP_RUN_SIMULATIONS = "true")
source("simulations/time_ar1/run_time.R")
```

The progress bar counts completed method fits. Development-time AR(1)
diagnostics are under `simulations/other_simus/ar1_method_comparisons/`.
