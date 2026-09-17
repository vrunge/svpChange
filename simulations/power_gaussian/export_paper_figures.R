## Export the Gaussian Section 5.1 figures to SVP_NEW_Figures2.

library(dplyr)
library(ggplot2)

source(file.path("simulations", "power_common.R"))

GAUSSIAN_FIGURE_ROOT <- normalizePath(
  file.path("..", "SVP_NEW_Figures2"), mustWork = FALSE
)

export_gaussian_paper_figures <- function(
    results = readRDS(file.path("simulations", "power_gaussian", "results.rds")),
    root = GAUSSIAN_FIGURE_ROOT) {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  scenario <- power_scenario_plot(1000L, 0.6, stats::rnorm)
  mse <- metric_plot(results, "MSE", "Signal MSE")
  distribution <- changepoint_distribution_plot(
    results, n = 1000L, selected_jump = 0.6
  )
  plots <- list(
    "1_signal_scenarios.pdf" = scenario,
    "1_f1_vs_jump.pdf" = metric_plot(results, "F1", "F1", TRUE),
    "1_precision_vs_jump.pdf" =
      metric_plot(results, "Precision", "Precision", TRUE),
    "1_recall_vs_jump.pdf" = metric_plot(results, "Recall", "Recall", TRUE),
    "1_correct_number_vs_jump.pdf" = metric_plot(
      results, "CorrectNumCP",
      "Probability of the correct number of changepoints", TRUE
    ),
    "1_numsegments_vs_jump.pdf" = metric_plot(
      results, "NumSegments", "Number of segments"
    ),
    "1_mse_vs_jump.pdf" = mse,
    "1_changepoint_distributions.pdf" = distribution
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

if (identical(tolower(Sys.getenv("SVP_EXPORT_GAUSSIAN_FIGURES")), "true")) {
  export_gaussian_paper_figures()
}
