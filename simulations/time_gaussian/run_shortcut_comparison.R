## Re-run the canonical Gaussian benchmark after the SVP shortcut changes.
## The checked-in gaussian_time_* files are read-only baseline data.

previous_run_flag <- Sys.getenv("SVP_RUN_SIMULATIONS", unset = NA_character_)
Sys.setenv(SVP_RUN_SIMULATIONS = "false")
source(file.path("simulations", "time_gaussian", "run_time.R"))
if (is.na(previous_run_flag)) {
  Sys.unsetenv("SVP_RUN_SIMULATIONS")
} else {
  Sys.setenv(SVP_RUN_SIMULATIONS = previous_run_flag)
}

SHORTCUT_PREFIX <- "gaussian_time_shortcuts"
KEY_COLUMNS <- c("experiment", "n", "k", "rep", "method")

pair_gaussian_timings <- function(baseline, shortcuts) {
  if (anyDuplicated(baseline[KEY_COLUMNS]) ||
      anyDuplicated(shortcuts[KEY_COLUMNS])) {
    stop("The timing key is not unique in one of the result tables.")
  }

  paired <- merge(
    baseline, shortcuts,
    by = KEY_COLUMNS,
    suffixes = c("_baseline", "_shortcuts"),
    all = TRUE,
    sort = TRUE
  )
  if (nrow(paired) != nrow(baseline) || nrow(paired) != nrow(shortcuts) ||
      anyNA(paired[c("time_baseline", "time_shortcuts")])) {
    stop("The baseline and shortcut benchmark keys do not agree.")
  }
  paired$detected_equal <-
    paired$detected_baseline == paired$detected_shortcuts

  # The checked-in benchmark predates the subtests API and the subsequent
  # multiscale-candidate corrections. Treat it as a historical comparison,
  # not as the shortcut regression oracle. PELT and both right-mode methods
  # must nevertheless reproduce their historical segment counts exactly.
  stable_methods <- paired$method != "SVP BIC multiscale"
  if (any(!paired$detected_equal[stable_methods])) {
    stop(
      "Detected segment counts changed for PELT or a right-mode method in ",
      sum(!paired$detected_equal[stable_methods]), " benchmark rows."
    )
  }

  paired$speedup <- ifelse(
    paired$time_baseline > 0 & paired$time_shortcuts > 0,
    paired$time_baseline / paired$time_shortcuts,
    NA_real_
  )
  paired
}

summarize_paired_timings <- function(paired) {
  groups <- split(
    paired,
    interaction(
      paired$experiment, paired$n, paired$k, paired$method,
      drop = TRUE, lex.order = TRUE
    )
  )
  dplyr::bind_rows(lapply(groups, function(group) {
    baseline_median <- stats::median(group$time_baseline)
    shortcut_median <- stats::median(group$time_shortcuts)
    data.frame(
      experiment = group$experiment[1L],
      n = group$n[1L],
      k = group$k[1L],
      method = group$method[1L],
      baseline_median = baseline_median,
      shortcut_median = shortcut_median,
      median_speedup = if (
        baseline_median > 0 && shortcut_median > 0
      ) baseline_median / shortcut_median else NA_real_,
      matching_segment_counts = sum(group$detected_equal),
      paired_rows = nrow(group)
    )
  }))
}

plot_shortcut_comparison <- function(paired, summary, root) {
  plot_root <- file.path(root, "plots")
  dir.create(plot_root, recursive = TRUE, showWarnings = FALSE)

  positive <- dplyr::filter(
    paired, time_baseline > 0, time_shortcuts > 0
  )
  before_after <- ggplot2::ggplot(
    positive,
    ggplot2::aes(time_baseline, time_shortcuts, colour = method)
  ) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey50") +
    ggplot2::geom_point(alpha = 0.25, size = 0.8) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::coord_equal() +
    ggplot2::labs(
      x = "Baseline elapsed time (s)",
      y = "Version 0.2.0 elapsed time (s)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top", legend.direction = "horizontal")

  speedup_n <- ggplot2::ggplot(
    dplyr::filter(
      summary, experiment == "vary_n", is.finite(median_speedup)
    ),
    ggplot2::aes(n, median_speedup, colour = method)
  ) +
    ggplot2::geom_hline(yintercept = 1, colour = "grey50") +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = "Sequence length", y = "Median baseline/new speedup") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top", legend.direction = "horizontal")

  speedup_k <- ggplot2::ggplot(
    dplyr::filter(
      summary, experiment == "vary_k", is.finite(median_speedup)
    ),
    ggplot2::aes(k, median_speedup, colour = method)
  ) +
    ggplot2::geom_hline(yintercept = 1, colour = "grey50") +
    ggplot2::geom_smooth(method = "loess", se = FALSE) +
    ggplot2::labs(x = "True number of changes", y = "Median baseline/new speedup") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top", legend.direction = "horizontal")

  ggplot2::ggsave(
    file.path(plot_root, paste0(SHORTCUT_PREFIX, "_before_after.pdf")),
    before_after, width = 7, height = 6
  )
  ggplot2::ggsave(
    file.path(plot_root, paste0(SHORTCUT_PREFIX, "_speedup_vs_n.pdf")),
    speedup_n, width = 8, height = 5
  )
  ggplot2::ggsave(
    file.path(plot_root, paste0(SHORTCUT_PREFIX, "_speedup_vs_k.pdf")),
    speedup_k, width = 8, height = 5
  )
  invisible(list(
    before_after = before_after,
    speedup_vs_n = speedup_n,
    speedup_vs_k = speedup_k
  ))
}

save_shortcut_comparison <- function(baseline, shortcuts) {
  paired <- pair_gaussian_timings(baseline, shortcuts)
  summary <- summarize_paired_timings(paired)
  utils::write.csv(
    paired,
    file.path(
      GAUSSIAN_TIME_ROOT,
      paste0(SHORTCUT_PREFIX, "_comparison.csv")
    ),
    row.names = FALSE
  )
  utils::write.csv(
    summary,
    file.path(
      GAUSSIAN_TIME_ROOT,
      paste0(SHORTCUT_PREFIX, "_summary.csv")
    ),
    row.names = FALSE
  )
  plot_shortcut_comparison(paired, summary, GAUSSIAN_TIME_ROOT)
  message(
    "Historical detected segment counts agree in ",
    sum(paired$detected_equal), " of ", nrow(paired), " paired fits; ",
    sum(!paired$detected_equal),
    " multiscale rows differ from the pre-subtests baseline."
  )
  invisible(list(results = shortcuts, comparison = paired, summary = summary))
}

run_shortcut_comparison <- function() {
  baseline_path <- file.path(
    GAUSSIAN_TIME_ROOT, "gaussian_time_results.csv"
  )
  if (!file.exists(baseline_path)) {
    stop("Missing checked-in baseline: ", baseline_path)
  }

  baseline <- utils::read.csv(baseline_path, stringsAsFactors = FALSE)
  shortcuts <- run_time_gaussian(
    reps_n = 100L, reps_k = 100L,
    reps_none_n = 100L, reps_none_k = 100L
  )
  save_time_outputs(shortcuts, GAUSSIAN_TIME_ROOT, SHORTCUT_PREFIX)
  save_shortcut_comparison(baseline, shortcuts)
}

if (identical(tolower(Sys.getenv("SVP_RUN_SIMULATIONS")), "true")) {
  run_shortcut_comparison()
}
