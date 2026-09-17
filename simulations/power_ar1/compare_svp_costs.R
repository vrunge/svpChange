## Threshold-matched comparison of the two SVP within-segment costs over the full
## Section 5.3 power grid. Run after run_ar1_calibration(); writes
## svp_cost_comparison.csv. Both variants use the calibrated SVP threshold, so
## the only thing that differs is how the partitions with the smallest number of
## valid segments are ranked.
suppressMessages(library(dplyr))
source(file.path("simulations", "power_ar1", "run_power.R"))
cal <- readRDS(file.path(AR1_ROOT, "ar1_calibration.rds"))
cst <- cal$selected$svp_constant

res <- run_power_grid(
  AR1_N, AR1_JUMP_SIZES, AR1_REPS,
  simulate_noise = function(size) simulate_ar1_noise(size, AR1_RHO),
  fit_methods = function(y, ...) list(
    ar1      = svp_ar1focus_boundaries(y, AR1_RHO, constant = cst, cost = "ar1"),
    gaussian = svp_ar1focus_boundaries(y, AR1_RHO, constant = cst, cost = "gaussian")),
  tolerance = AR1_TOLERANCE, workers = power_default_workers(), seed = 123L)
res <- complete_localization_error(res)

loc <- res |> filter(CorrectNumCP == 1, pattern != "none",
                     is.finite(LocalizationError)) |>
  group_by(pattern, cost = algorithm) |>
  summarise(LocalizationError = mean(LocalizationError), .groups = "drop")
tab <- res |> group_by(pattern, cost = algorithm) |>
  summarise(F1 = mean(F1), Precision = mean(Precision), Recall = mean(Recall),
            CorrectNumCP = mean(CorrectNumCP), MSE = mean(MSE),
            NumSegments = mean(NumSegments), .groups = "drop") |>
  left_join(loc, by = c("pattern", "cost"))
tab$constant <- cst
tab$reps <- AR1_REPS
tab$null_target_f1 <- cal$selected$target_f1
tab$decafs_version <- cal$decafs_version
write.csv(tab, file.path(AR1_ROOT, "svp_cost_comparison.csv"), row.names = FALSE)
print(as.data.frame(tab[order(tab$pattern, tab$cost), ]), digits = 4)
cat("\nsegment counts identical for every single series:",
    isTRUE(all.equal(res$NumSegments[res$algorithm == "ar1"],
                     res$NumSegments[res$algorithm == "gaussian"])), "\n")
