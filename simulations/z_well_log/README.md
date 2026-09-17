# Well-log SVP calibration review

This directory is an isolated review artifact. It does not replace or modify
the existing well-log application or the manuscript.

From the package root, with the current package source installed:

```r
source("simulations/z_well_log/run_review.R")
review <- run_well_log_review(workers = 8L)
rmarkdown::render(
  "simulations/z_well_log/review_well_log.Rmd",
  output_file = "review_well_log.pdf",
  output_dir = "simulations/z_well_log",
  knit_root_dir = "simulations/z_well_log"
)
```

The recommended fit is Median--Mood SVP with `subtests = "right"` and a
threshold calibrated to 0.99 one-segment recovery under a residual moving-block
bootstrap. The calibration uses the maximum validity statistic over every
growing prefix, which is the exact one-segment condition for right pruning.
RFPOP is used only as an external comparison segmentation.

The `both` mode is included as a sensitivity analysis. It is not labelled as
calibrated by the right-mode threshold because any invalid interior subsegment
can prune the one-segment candidate. `both_null_stress.csv` records a direct
null check of that distinction.
