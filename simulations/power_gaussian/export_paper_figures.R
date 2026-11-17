## Export the Gaussian Section 5.1 figures to simulations/all_plots.

library(dplyr)
library(ggplot2)

source(file.path("simulations", "power_common.R"))
source(file.path("simulations", "gaussian_common.R"))

GAUSSIAN_FIGURE_ROOT <- normalizePath(
  file.path("simulations", "all_plots"), mustWork = FALSE
)

export_gaussian_paper_figures <- function(
    results = readRDS(file.path("simulations", "power_gaussian", "results.rds")),
    root = GAUSSIAN_FIGURE_ROOT, include_scenario = TRUE) {
  order <- validate_gaussian_power_results(results)
  results <- set_algorithm_order(results, order)
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  mse <- metric_plot(results, "MSE", "Signal MSE")
  distribution <- changepoint_distribution_plot(
    results, n = 1000L, selected_jump = 0.6
  )
  plots <- list(
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
  if (include_scenario) {
    plots[["1_signal_scenarios.pdf"]] <- power_scenario_plot(
      1000L, 0.6, stats::rnorm
    )
  }
  for (file in names(plots)) {
    ggplot2::ggsave(
      file.path(root, file), plots[[file]],
      width = if (grepl("distributions", file)) 15 else
        if (grepl("signal_scenarios", file)) 10 else 12,
      height = if (grepl("distributions", file)) 12 else 7
    )
  }
  four <- readRDS(file.path("simulations", "power_gaussian",
                                "f1_four_methods_results.rds"))
  four_labels <- levels(four$algorithm)
  four_colours <- stats::setNames(
    scales::hue_pal()(length(four_labels)), four_labels
  )
  four_plot <- metric_plot(four, "F1", "F1", TRUE) +
    ggplot2::scale_colour_manual(values = four_colours) +
    ggplot2::scale_fill_manual(values = four_colours) +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 3L, byrow = TRUE)) +
    ggplot2::theme(legend.text = ggplot2::element_text(size = 7.5),
                   legend.title = ggplot2::element_text(size = 9),
                   legend.key.width = grid::unit(0.5, "lines"),
                   legend.spacing.x = grid::unit(0.25, "lines"))
  ggplot2::ggsave(file.path(root, "f1_four_methods.pdf"), four_plot,
                  width = 14, height = 8, limitsize = FALSE)

  ## Subset figure with only the three original methods (PELT, SVP both,
  ## SVP SMUCE), which gaussian_four_method_labels() always places first.
  three_labels <- four_labels[1:3]
  three <- dplyr::filter(four, algorithm %in% three_labels)
  three$algorithm <- factor(as.character(three$algorithm), levels = three_labels)
  ## Standard ggplot2 default palette for exactly these three groups (not
  ## a slice of the eleven-method palette).
  three_plot <- metric_plot(three, "F1", "F1", TRUE) +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 1L, byrow = TRUE))
  ggplot2::ggsave(file.path(root, "f1_three_methods.pdf"), three_plot,
                  width = 12, height = 7)

  ## Subset figure with the three original methods plus the three
  ## best-performing "both"-family variants (by mean F1 across all
  ## scenarios and jump sizes, and separately in the harder jump<=1.0
  ## regime: both_rmin, both_correction, both_mosum, in that order both
  ## times -- see multiscale_methods_tested.pdf).
  six_labels <- c(three_labels,
                  four_labels[grepl("^SVP both_rmin ", four_labels)],
                  four_labels[grepl("^SVP both_correction ", four_labels)],
                  four_labels[grepl("^SVP both_mosum \\(", four_labels)])
  six <- dplyr::filter(four, algorithm %in% six_labels)
  six$algorithm <- factor(as.character(six$algorithm), levels = six_labels)
  ## Standard ggplot2 six-colour default, but reordered so the three
  ## methods shared with f1_three_methods.pdf keep EXACTLY the same colour
  ## there (hue_pal(3) is always positions 1,3,5 of hue_pal(6), since 6 is
  ## a multiple of 3): PELT/SVP both/SVP SMUCE match f1_three_methods.pdf
  ## pixel-for-pixel; both_rmin/both_correction/both_mosum take the three
  ## remaining standard hues.
  six_hues <- scales::hue_pal()(6)[c(1, 3, 5, 2, 4, 6)]
  six_colours <- stats::setNames(six_hues, six_labels)
  six_plot <- metric_plot(six, "F1", "F1", TRUE) +
    ggplot2::scale_colour_manual(values = six_colours) +
    ggplot2::scale_fill_manual(values = six_colours) +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2L, byrow = TRUE))
  ggplot2::ggsave(file.path(root, "f1_six_methods.pdf"), six_plot,
                  width = 13, height = 7.5)
  invisible(plots)
}

if (identical(tolower(Sys.getenv("SVP_EXPORT_GAUSSIAN_FIGURES")), "true")) {
  export_gaussian_paper_figures()
}
