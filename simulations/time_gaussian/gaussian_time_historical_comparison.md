# Gaussian timing audit

- The historical Experiment 1 calls `SVP()` without an explicit mode;
  its meaning depends on the package version installed when it was run.
- The checked-in paper timing script currently assigns the same `both`,
  `c = 1.5` call to two differently named curves. Its saved CSV does not
  show identical timings for those curves, so script and result provenance
  are not internally reproducible as checked in.
- The refreshed script specifies every mode and constant, warms each method,
  randomizes method order, shares data, and batches sub-resolution calls.
- Absolute old/new time ratios are not treated as hardware-independent.
- Audited paper script: `../svpChange_old/simulations_paper/time_gaussian/run_time.R`.
