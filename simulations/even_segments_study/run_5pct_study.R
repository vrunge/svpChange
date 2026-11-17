## Same even-segment up--down study with the 5% calibration level.
##
## Uses the same simulated series as the 1% study (same deterministic seeds)
## so that the two calibration levels are paired. Constants are read from
## the alpha study: PELT 1.725 and SVP both 3.604 at the 5% target.
##
## Run from the package root:
##   SVP_RUN_EVEN_5PCT=true Rscript simulations/even_segments_study/run_5pct_study.R

source(file.path("simulations", "power_common.R"))
source(file.path("simulations", "even_segments_study", "run_study.R"))

library(ggplot2)

read_5pct_constants <- function() {
  calibration <- utils::read.csv(file.path(
    "simulations", "even_segments_study", "calibration_5pct.csv"
  ))
  pelt <- calibration$constant[
    calibration$alpha == 0.05 & calibration$method == "PELT calibrated"
  ]
  svp <- calibration$constant[
    calibration$alpha == 0.05 & grepl("SVP both", calibration$method)
  ]
  stopifnot(length(pelt) == 1L, length(svp) == 1L)
  list(pelt = pelt, svp = svp)
}

plot_segments_curves <- function(results, path) {
  results$pelt_segments <- results$pelt_cp_count + 1
  results$svp_segments <- results$svp_cp_count + 1
  pelt_label <- results$pelt_label[[1L]]
  svp_label <- results$svp_label[[1L]]
  curves <- rbind(
    data.frame(segments = results$segments, jump = results$jump,
               algorithm = pelt_label, segments_hat = results$pelt_segments),
    data.frame(segments = results$segments, jump = results$jump,
               algorithm = svp_label, segments_hat = results$svp_segments)
  )
  curves$algorithm <- factor(curves$algorithm, c(pelt_label, svp_label))
  summary <- aggregate(segments_hat ~ segments + jump + algorithm, curves, mean)
  se <- aggregate(segments_hat ~ segments + jump + algorithm, curves,
                  function(v) stats::sd(v) / sqrt(length(v)))
  summary$se <- se$segments_hat
  summary$lower <- pmax(0, summary$segments_hat - 1.96 * summary$se)
  summary$upper <- summary$segments_hat + 1.96 * summary$se
  summary$panel <- factor(sprintf("%d segments", summary$segments),
                          sprintf("%d segments", EVEN_SEGMENTS))
  truth <- data.frame(
    panel = factor(sprintf("%d segments", EVEN_SEGMENTS),
                   sprintf("%d segments", EVEN_SEGMENTS)),
    true_segments = EVEN_SEGMENTS
  )
  plot <- ggplot2::ggplot(
    summary,
    ggplot2::aes(jump, segments_hat, colour = algorithm, fill = algorithm)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.12, colour = NA
    ) +
    ggplot2::geom_hline(
      data = truth, ggplot2::aes(yintercept = true_segments),
      linetype = 3, colour = "grey30", inherit.aes = FALSE
    ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 0.9) +
    ggplot2::facet_wrap(~panel, nrow = 2L, scales = "free_y") +
    ggplot2::labs(
      x = "Jump size", y = "Mean number of estimated segments",
      colour = "Method", fill = "Method"
    ) +
    paper_theme()
  ggplot2::ggsave(path, plot, width = 16, height = 7)
  invisible(plot)
}

run_5pct_study <- function(reps = EVEN_REPS, workers = power_default_workers(),
                           root = EVEN_ROOT) {
  constants <- read_5pct_constants()
  results <- dplyr::bind_rows(lapply(seq_along(EVEN_SEGMENTS), function(i) {
    cat("segments", EVEN_SEGMENTS[[i]], "\n")
    even_batch(i, reps, constants$pelt, constants$svp, workers)
  }))
  results$pelt_label <- sprintf("PELT calibrated, c = %.3f", constants$pelt)
  results$svp_label <- sprintf("SVP both calibrated, c = %.3f", constants$svp)
  saveRDS(results, file.path(root, "results_5pct.rds"), compress = FALSE)
  plot_f1_curves(results, file.path("simulations", "all_plots",
                                    "6_f1_vs_segments_5pct.pdf"))
  plot_segments_curves(results, file.path("simulations", "all_plots",
                                          "6_segments_vs_jump_5pct.pdf"))
  results
}

if (identical(tolower(Sys.getenv("SVP_RUN_EVEN_5PCT")), "true")) {
  results <- run_5pct_study()
  means <- aggregate(cbind(pelt_f1, svp_f1) ~ segments, results, mean)
  cat("\nMean F1 by segment count (5% calibration):\n")
  print(means, row.names = FALSE, digits = 3)
}
