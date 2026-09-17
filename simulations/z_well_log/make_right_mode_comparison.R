## Compare the two proposed right-pruned SVP configurations.

WELL_LOG_ROOT <- file.path("simulations", "z_well_log")
PLOT_FILE <- file.path(
  WELL_LOG_ROOT, "plots", "06_wilcoxon_mood_right_comparison.pdf"
)
CSV_FILE <- file.path(
  WELL_LOG_ROOT, "wilcoxon_mood_right_comparison.csv"
)

source(file.path(WELL_LOG_ROOT, "run_review.R"))
require_well_log_packages()
review <- readRDS(file.path(WELL_LOG_ROOT, "well_log_review.rds"))
y <- as.numeric(DeCAFS::oilWell)
n <- length(y)
reference <- review$reference

wilcoxon_K <- 10L
wilcoxon_gamma <- 1.5 * sqrt((n / wilcoxon_K)^3 / 12)
mood_gamma <- review$metadata$selected_gamma

fits <- list(
  reference,
  fit_svp(y, wilcoxon_gamma, "WilcoxonCost", "right"),
  fit_svp(y, mood_gamma, "MedianMoodCost", "right")
)
names(fits) <- c(
  "RFPOP benchmark\nbiweight loss, lambda = 70, lthreshold = 2",
  paste0(
    "Wilcoxon SVP / right\nK = ", wilcoxon_K,
    ", gamma = ", formatC(wilcoxon_gamma, format = "f", digits = 1)
  ),
  paste0(
    "Median Mood SVP / right\ncalibrated gamma = ",
    formatC(mood_gamma, format = "f", digits = 2)
  )
)
method_levels <- names(fits)

summary_rows <- lapply(names(fits), function(method) {
  changepoints <- internal_boundaries(fits[[method]], n)
  score <- match_changepoints(
    reference, changepoints, WELL_LOG_MATCH_TOLERANCE
  )
  data.frame(
    method = method,
    n_changepoints = length(changepoints),
    matched = unname(score[["tp"]]),
    precision = unname(score[["precision"]]),
    recall = unname(score[["recall"]]),
    F1 = unname(score[["F1"]]),
    changepoints = paste(changepoints, collapse = ","),
    stringsAsFactors = FALSE
  )
})
summary_table <- do.call(rbind, summary_rows)
utils::write.csv(summary_table, CSV_FILE, row.names = FALSE)

plot_data <- do.call(rbind, lapply(names(fits), function(method) {
  changepoints <- internal_boundaries(fits[[method]], n)
  fitted <- segment_median_fit(y, changepoints)$fitted
  data.frame(
    position = seq_len(n),
    observation = y,
    fitted = fitted,
    method = method,
    stringsAsFactors = FALSE
  )
}))
plot_data$method <- factor(plot_data$method, levels = method_levels)

plot <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(position, observation)
) +
  ggplot2::geom_point(size = 0.18, alpha = 0.45) +
  ggplot2::geom_line(
    ggplot2::aes(y = fitted), linewidth = 0.7, colour = "#AA3377"
  ) +
  ggplot2::geom_vline(
    data = do.call(rbind, lapply(names(fits), function(method) {
      data.frame(
        position = internal_boundaries(fits[[method]], n),
        method = factor(method, levels = method_levels),
        stringsAsFactors = FALSE
      )
    })),
    ggplot2::aes(xintercept = position),
    linetype = 3, linewidth = 0.3, colour = "grey35"
  ) +
  ggplot2::facet_wrap(~method, ncol = 1) +
  ggplot2::scale_x_continuous(
    breaks = seq(0, n, by = 1000), expand = ggplot2::expansion(mult = c(0.01, 0.01))
  ) +
  ggplot2::labs(
    x = "Sequence position", y = "Log measurement",
    caption = paste0(
      "Vertical lines are detected changepoints; purple curves are " ,
      "within-segment medians. Matching tolerance: ",
      WELL_LOG_MATCH_TOLERANCE, " observations."
    )
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(face = "bold", hjust = 0),
    axis.text.x = ggplot2::element_text(size = 9),
    plot.caption = ggplot2::element_text(hjust = 0, size = 8),
    panel.spacing = grid::unit(0.8, "lines")
  )

ggplot2::ggsave(PLOT_FILE, plot, width = 10, height = 8.5)
print(summary_table)
