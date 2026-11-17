## Export corrected Section 5.3 figures to simulations/all_plots.

library(ggplot2)

source(file.path("simulations", "power_common.R"))
source(file.path("simulations", "power_ar1", "calibrate_ar1.R"))

AR1_FIGURE_ROOT <- normalizePath(
  file.path("simulations", "all_plots"), mustWork = FALSE
)

export_ar1_paper_figures <- function(
    results = readRDS(file.path("simulations", "power_ar1", "results.rds")),
    calibration = readRDS(file.path(
      "simulations", "power_ar1", "ar1_calibration.rds"
    )),
    root = AR1_FIGURE_ROOT,
    selected_jump = 1.4, include_scenario = TRUE) {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  results <- complete_localization_error(results)
  results <- set_algorithm_order(results, ar1_method_labels(calibration))
  plots <- list(
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
  if (include_scenario) {
    plots[["5_signal_scenarios_ar1.pdf"]] <- power_scenario_plot(
      AR1_N, selected_jump,
      function(n) simulate_ar1_noise(n, AR1_RHO)
    )
  }
  ## The five AR(1) method labels are long; a single-row legend overflows
  ## the plot width, so wrap it onto two rows.
  legend_two_rows <- ggplot2::guides(
    colour = ggplot2::guide_legend(nrow = 2L, byrow = TRUE),
    fill = ggplot2::guide_legend(nrow = 2L, byrow = TRUE)
  )
  for (file in names(plots)) {
    if (!grepl("distributions|signal_scenarios", file)) {
      plots[[file]] <- plots[[file]] + legend_two_rows
    }
  }
  for (file in names(plots)) {
    ggplot2::ggsave(
      file.path(root, file), plots[[file]],
      width = if (grepl("distributions", file)) 15 else
        if (grepl("signal_scenarios", file)) 10 else 12,
      height = if (grepl("distributions", file)) 12 else 7
    )
  }
  invisible(plots)
}

if (identical(tolower(Sys.getenv("SVP_EXPORT_AR1_FIGURES")), "true")) {
  export_ar1_paper_figures()
}
