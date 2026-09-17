# Corrected Student-t(2) Section 5.2 power study

Run from the `svpChange` package root:

```r
source("simulations/power_robust/calibrate_robust.R")
calibration <- run_robust_calibration(workers = 8L)
source("simulations/power_robust/run_power.R")
run_and_save_robust(workers = 8L)
```

The default calibration uses 20,000 independent Student-t(2) null series and
10,000 independent validation series. Seeds, candidate values, Wilson
intervals, the selected parameters, and both PELT benchmarks are stored in
`true_true_calibration.csv` and `robust_calibration.rds`.

The RFPOP multiplier is the smallest candidate whose corrected null F1 is at
least 0.99. Wilcoxon and Median-Mood candidates are then selected by closest
corrected null F1 to that RFPOP target. The four headline algorithms are
fixed Gaussian PELT with penalty `2*log(n)`, calibrated RFPOP with biweight
loss and `lthreshold = 3`, `SVP MedianMood / right`, and `SVP Wilcoxon /
right`. There are no `both` or `multiscale` SVP calls in this study.

The Wilcoxon threshold is

```text
gamma_W(c, K) = c * sqrt((n/K)^3 / 12)
```

and the Mood threshold uses the Pearson chi-square law with
`alpha_s = 1 - (1 - alpha)^(1/(n/K - 1))`. Here `K` is passed to the fitting
code as `oracle_segments`: it is the known data-generating segment count in
this power experiment, not an available input in practical use.

Metrics use the shared corrected framework: the mandatory endpoint `n` is
removed, tolerance matches are one-to-one, and F1, precision, recall,
correct-number probability, MSE, segment count, and localization error are
reported. The Gaussian-PELT null F1 is an external benchmark from the
corrected Gaussian evaluation; the heavy-tail PELT fit is intentionally
reported separately because its misspecified Gaussian SSE fit has almost
universal false positives under t(2).
