## Export the corrected robust Section 5.2 figures to simulations/all_plots.

library(svpChange)
library(dplyr)
library(ggplot2)

source(file.path("simulations", "power_common.R"))
source(file.path("simulations", "power_robust", "calibrate_robust.R"))

ROBUST_ROOT <- file.path("simulations", "power_robust")
FIGURE_ROOT <- file.path("simulations", "all_plots")

export_robust_paper_figures <- function(
    results = readRDS(file.path(ROBUST_ROOT, "results.rds")),
    calibration = readRDS(file.path(ROBUST_ROOT, "robust_calibration.rds")),
    root = FIGURE_ROOT, selected_jump = 0.6,
    include_scenario = TRUE) {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  results <- set_algorithm_order(results, robust_method_labels(calibration))
  results <- complete_localization_error(results)
  precision <- metric_plot(results, "Precision", "Precision", TRUE)
  recall <- metric_plot(results, "Recall", "Recall", TRUE)
  distribution <- changepoint_distribution_plot(
    results, 1000L, selected_jump
  )
  plots <- list(
    "3_f1_vs_jump.pdf" = metric_plot(results, "F1", "F1", TRUE),
    "3_precision_vs_jump.pdf" = precision,
    "3_recall_vs_jump.pdf" = recall,
    "3_correct_number_vs_jump.pdf" = metric_plot(
      results, "CorrectNumCP",
      "Probability of the correct number of changepoints", TRUE
    ),
    "3_numsegments_vs_jump.pdf" = metric_plot(
      results, "NumSegments", "Number of segments"
    ),
    "3_mse_vs_jump.pdf" = metric_plot(
      results, "MSE", "Signal MSE"
    ),
    "3_localizationerror_vs_jump.pdf" = localization_error_plot(results),
    "3_changepoint_distributions_robust.pdf" = distribution
  )
  if (include_scenario) {
    plots[["3_signal_scenarios_robust.pdf"]] <- power_scenario_plot(
      1000L, selected_jump,
      function(n) ts_generator(
        chpts = n, parameters = 0, type = "student", df = 2, scale = 1
      ),
      y_limits = c(-10, 10)
    )
  }
  sizes <- lapply(names(plots), function(file) {
    if (grepl("distributions", file)) c(15, 12) else
      if (grepl("localizationerror", file)) c(12, 7) else c(10, 7)
  })
  for (i in seq_along(plots)) {
    ggsave(
      file.path(root, names(plots)[i]), plots[[i]],
      width = sizes[[i]][1L], height = sizes[[i]][2L]
    )
  }
  writeLines(c(
    "Corrected Student-t(2) robust Section 5.2 paper figures.",
    paste0("Selected RFPOP c = ", calibration$selected$rfpop_constant),
    paste0("Selected Wilcoxon c = ", calibration$selected$wilcoxon_constant),
    paste0("Selected MedianMood alpha = ", calibration$selected$mood_alpha),
    "All SVP calls use subtests = right; K is the oracle data-generating segment count.",
    paste0("Calibration target F1 = ", calibration$selected$rfpop_target_f1),
    paste0("Gaussian-PELT benchmark F1 = ",
           calibration$gaussian_pelt_benchmark$gaussian_pelt_null_f1)
  ), file.path(root, "robust_figure_manifest.txt"))
  invisible(plots)
}

if (identical(tolower(Sys.getenv("SVP_EXPORT_ROBUST_FIGURES")), "true")) {
  export_robust_paper_figures()
}
