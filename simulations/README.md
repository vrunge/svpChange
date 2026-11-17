# Figures used by the paper

`all_plots/` is the single destination for the paper PDFs. Run scripts from
the package root. Each study's calibration and power or timing driver remains
here only when it contributes to a PDF in `all_plots/`.

| PDFs in `all_plots/` | Producer | Saved input, when available |
| --- | --- | --- |
| `1_*`, `f1_four_methods.pdf`, `updown_four_jump_sizes.pdf` | `power_gaussian/export_paper_figures.R`, `power_gaussian/plot_updown_jump_sizes.R` | `power_gaussian/results.rds`, `power_gaussian/f1_four_methods_results.rds` |
| `2_gaussian_time.pdf` | `time_gaussian/run_time.R` | `time_gaussian/gaussian_time_results.rds` |
| `2_robust_time.pdf` | `time_robust/run_time.R` | `time_robust/robust_time_results.rds` |
| `2_ar1_time.pdf` | `time_ar1/run_time.R` | `time_ar1/ar1_time_results.rds` |
| `3_*` | `power_robust/export_paper_figures.R` | Calibration in `power_robust/robust_calibration.rds` |
| `4_segmentation_logdata.pdf` | `well_log/run.R` | `DeCAFS::oilWell` |
| `5_*` | `power_ar1/export_paper_figures.R` | Calibration in `power_ar1/ar1_calibration.rds` |
| `6_*` | `even_segments_study/run_study.R`, `run_5pct_study.R`, `make_segments_figure.R`, `make_example_8seg.R` | `even_segments_study/results.rds`, `results_5pct.rds` |

The Gaussian mean-change statistic is the **full maximum RSS reduction**.
Known-variance AR(1) validity uses the full innovation RSS reduction divided
by `sigma2`. The matching SVP gamma values in saved data and figure labels are
on those scales. PELT penalties and the separate multiscale statistics retain
their own definitions.

All figures in `all_plots/` are regenerated from fresh runs against the
current package: the Gaussian, AR(1), and robust power studies were rerun
end to end (`power_gaussian/run_power.R`, `power_ar1/run_power.R`,
`power_robust/run_power.R`), and the timing, well-log, and even-segment
panels were redrawn from their saved result objects.
