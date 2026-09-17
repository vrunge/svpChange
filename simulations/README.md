# Reproducible simulation studies

Run every script from the `svpChange` package root. The six primary studies
are deliberately arranged as three power studies and three runtime studies:

- `power_gaussian/`
- `power_robust/`
- `power_ar1/`
- `time_gaussian/`
- `time_robust/`
- `time_ar1/`

Each primary folder has one entry point (`run_power.R` or `run_time.R`), one
README, and its generated results. Supporting experiments and historical
scripts live under `other_simus/`.

The paper-figure refresh is driven by
`other_simus/maintenance/rebuild_svp_new_figures2.R`. It writes the 33
updated PDFs (and current result/report files) to the sibling folder
`../SVP_NEW_Figures2`, leaving the historical `../SVP_NEW_Figures` untouched.
The default quick mode consumes the maintained v0.2.0 result objects; set
`SVP_REBUILD_FULL=true` together with `SVP_REBUILD_FIGURES2=true` to run the
approved 1000-null/200-signal calibration grids and regenerate all power
studies.

Shared code is defined once:

- `power_common.R` contains power-study simulation, metrics, plotting, and
  saving infrastructure.
- `time_common.R` contains timing grids, progress reporting, signal helpers,
  plotting, and saving infrastructure.

Method names stored in result files are the same names used in plot legends.
Every stochastic task receives a deterministic seed based on the study seed
and its design-row index. Consequently, sequential and fork-parallel power
runs generate the same simulated datasets.

## Running a study

Set the opt-in environment variable and source one entry point:

```r
Sys.setenv(SVP_RUN_SIMULATIONS = "true")
source("simulations/power_gaussian/run_power.R")
source("simulations/time_gaussian/run_time.R")
```

For a quick or sequential run, source the script with simulations disabled and
call its `run_*()` function with smaller arguments. Power studies accept
`workers = 1L`; timing studies display a progress bar.

Each `run_power.R` writes `results.rds`, `results.csv`, and all seven PDFs in
its own `plots/` directory. Each `run_time.R` writes its RDS/CSV result pair
and its PDFs in its own `plots/` directory. Timing values equal to zero
at the operating-system timer resolution remain in the result files but are
excluded from the two logarithmic-time plots.

| Study folders | Entry point | PDFs produced |
|---|---|---|
| `power_gaussian`, `power_robust`, `power_ar1` | `run_power.R` | `01_scenarios.pdf` through `07_loc_error.pdf` |
| `time_gaussian` | `run_time.R` | `gaussian_time_time_pair.pdf` plus the two individual panels |
| `time_robust` | `run_time.R` | `robust_time_time_pair.pdf` plus the two individual panels |
| `time_ar1` | `run_time.R` | `ar1_time_time_pair.pdf` plus the two individual panels |

The AR(1) power study requires the optional `DeCAFS` package (version 3.3.6 or
later for the corrected Section 5.3 study). Robust power
comparisons additionally require `robseg`. The corrected heavy-tail study
uses the dedicated calibration driver in `power_robust/calibrate_robust.R`;
its four-method output, paired audit, and Section 5.2 replacement memo are
documented in `power_robust/README.md` and `z_robust/`.
The corrected AR(1) study uses `power_ar1/calibrate_ar1.R`; its four-method
output, calibration table, review, and Section 5.3 replacement memo are
documented in `power_ar1/README.md` and `z_AR1/`.
