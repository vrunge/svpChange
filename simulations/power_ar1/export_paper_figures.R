## Export corrected Section 5.3 figures to SVP_NEW_Figures2.

library(ggplot2)

source(file.path("simulations", "power_common.R"))
source(file.path("simulations", "power_ar1", "calibrate_ar1.R"))

AR1_FIGURE_ROOT <- normalizePath(
  file.path("..", "SVP_NEW_Figures2"), mustWork = FALSE
)

export_ar1_paper_figures <- function(
    results = readRDS(file.path("simulations", "power_ar1", "results.rds")),
    calibration = readRDS(file.path(
      "simulations", "power_ar1", "ar1_calibration.rds"
    )),
    root = AR1_FIGURE_ROOT,
    selected_jump = 1.4) {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  results <- complete_localization_error(results)
  results <- set_algorithm_order(results, ar1_method_labels(calibration))
  plots <- list(
    "5_signal_scenarios_ar1.pdf" = power_scenario_plot(
      AR1_N, selected_jump,
      function(n) simulate_ar1_noise(n, AR1_RHO)
    ),
    "5_f1_vs_jump.pdf" = metric_plot(results, "F1", "F1", TRUE),
    "5_precision_vs_jump.pdf" = metric_plot(
      results, "Precision", "Precision", TRUE
    ),
    "5_recall_vs_jump.pdf" = metric_plot(results, "Recall", "Recall", TRUE),
    "5_correct_number_vs_jump.pdf" = metric_plot(
      results, "CorrectNumCP",
      "Probability of the correct number of changepoints", TRUE
    ),
    "5_mse_vs_jump.pdf" = metric_plot(results, "MSE", "Signal MSE"),
    "5_localizationerror_vs_jump.pdf" = localization_error_plot(results),
    "5_changepoint_distributions_ar1.pdf" = changepoint_distribution_plot(
      results, AR1_N, selected_jump
    )
  )
  for (file in names(plots)) {
    ggplot2::ggsave(
      file.path(root, file), plots[[file]],
      width = if (grepl("distributions", file)) 15 else 10,
      height = if (grepl("distributions", file)) 10 else 5
    )
  }
  invisible(plots)
}

if (identical(tolower(Sys.getenv("SVP_EXPORT_AR1_FIGURES")), "true")) {
  export_ar1_paper_figures()
}
