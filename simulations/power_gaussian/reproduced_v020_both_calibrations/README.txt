Current-source Gaussian extension using subtests = "both".

The independent null calibration uses n = 1000, 100 replications, seed = 700000,
and gamma = c log(n), with c from 1.4 to 2.4 by 0.1.  It gives null F1 = 0.973
at c = 1.8, 0.990 at c = 1.9, and 1.000 at c = 2.0.  The full Section 5.1
grid (seed = 123) evaluates c = 1.9 and c = 2.0; see results.csv and summary.csv.

The reproducible driver is run_both_calibrated.R.  The historical 1.8 multiplier
is retained in the canonical Section 5.1 reproduction for direct comparison.
The accompanying f1.pdf and null_calibration.pdf summarize the two comparisons.
