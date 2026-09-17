# Robust runtime study

`run_time.R` benchmarks the robust power-study competitors under Student-t(2)
noise. Method labels follow the power-study calibration:

- `PELT`
- `RFPOP`
- `SVP MedianMood / right / alpha = ...`
- `SVP Wilcoxon / right / c = ...`

The SVP constants are read from `power_robust/robust_calibration.rds`, and the
Mood and Wilcoxon thresholds use the simulated oracle number of segments.
This oracle-K convention is used only to measure the computational behavior
in the power-study design.

The sequence-length grid contains 16 values from 256 through 4,000; the
fixed-length experiment uses `n = 2000` and `k = 0, 2, ..., 40`. There are 64 replicates
per method and experiment. Fast fits are adaptively
batched for approximately 0.05 seconds before dividing by the batch size.
Methods are evaluated in a fixed order for reproducibility.

The preferred paper output is `robust_time_time_pair.pdf`: a side-by-side
two-panel figure with one shared legend above the panels. The individual
outputs remain `robust_time_time_vs_n.pdf` and
`robust_time_time_vs_detected.pdf` for compatibility. The first plot uses the median elapsed time
under the null; the second uses the median at each exact detected segment count
from the fixed-`n` experiment, with bootstrap 95% intervals. The old
`time_vs_k` file, if present, is historical and is not regenerated.

Run from the package root:

```r
Sys.setenv(SVP_RUN_SIMULATIONS = "true")
source("simulations/time_robust/run_time.R")
```

The separate C++ validity-update benchmark is under
`simulations/other_simus/robust_validity_benchmark/`.
