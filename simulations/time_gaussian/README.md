# Gaussian runtime study

`run_time.R` compares PELT and the three Gaussian SVP modes calibrated to the
same null false-positive probability:

- PELT / none / penalty = 2 log(n)
- SVP Gaussian / none / c = 0.968
- SVP Gaussian / right / c = 1.642
- SVP Gaussian / both / c = 1.987

The actual constants are always read from
`../power_gaussian/gaussian_mode_calibration.csv`. The benchmark varies the
sequence length under no change and the number of true changes at fixed
`n = 10000`. It uses `ts_generator()`, presents the same generated series to
all methods, randomizes their execution order, performs untimed warmups, and
adaptively batches fast calls to reduce timer-resolution error.

The one-hour paper run uses 12 replicates for PELT, SVP `right`, and SVP
`both`, and 48 replicates for SVP `none`. The `none` sample is larger because
its detected-segment curve is less stable. The sequence-length grid contains
24 values from 256 through 16,384; the fixed-length experiment uses
`n = 10000` and `k = 0:100`.

Run from the package root:

```r
Sys.setenv(SVP_RUN_SIMULATIONS = "false")
source("simulations/time_gaussian/run_time.R")
run_and_save_time_gaussian()
```

Raw per-call timings, median summaries, and bootstrap intervals are saved in
this folder. The detected-segment plot uses only the fixed-`n` experiment and
connects the median at each observed segment count; it does not use LOESS.
`gaussian_time_slopes.csv` reports descriptive log--log slopes for the complete
size grid and its upper half; they are empirical summaries, not worst-case
complexity proofs. The two paper-ready panels are copied to
`../SVP_NEW_Figures2/`.

Earlier implementations and timing results are retained under
`simulations/other_simus/time_gaussian_legacy/` for audit purposes.
