SMUCE/SVP COMPARISON

compare_smuce_svp.R compares `stepR::stepFit()` with the R and C++ constrained
Gaussian SVP implementations. It uses random heterogeneous Gaussian
piecewise-constant signals and compares the returned partitions directly.

Install the dependency first:

  install.packages("stepR")

Then run from the svpChange package root:

  source("simulations/smuce_comparison/compare_smuce_svp.R")
  answer <- compare_smuce_svp()
  answer$smuce
  answer$svp_smuce

The random seed is unset by default, so every run generates a new signal.
`run_smuce_replicates()` repeats the comparison and reports exact equality.

`time_smuce_svp.R` benchmarks `stepR::stepFit()` and
`svp_smuce_cpp()` for increasing sample sizes. The C++ routine updates
SMUCE envelopes incrementally in O(n^2) time. Both calls explicitly use all
subintervals, including for n > 1000. Timings depend on the signal and
hardware.

The comparison function returns:

  $smuce  the stepFit segment-end indices
  $svp_smuce  the constrained SVP segment-end indices
  $svp_smuce_cpp  the C++ constrained SVP segment-end indices

The concise derivation and exact argument mapping are in
simulations/smuce_comparison/SMUCE_SVP_equality.Rmd and its PDF.
