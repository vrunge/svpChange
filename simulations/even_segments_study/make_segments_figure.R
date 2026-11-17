## Ten-panel figure: mean number of estimated segments versus jump size,
## one panel per true segment count (2, 4, ..., 20), PELT versus SVP both.
## The dashed horizontal line in each panel is the true number of segments.
##
## Run from the package root:
##   Rscript simulations/even_segments_study/make_segments_figure.R

source(file.path("simulations", "power_common.R"))
source(file.path("simulations", "even_segments_study", "run_study.R"))

library(ggplot2)

results <- readRDS(file.path(EVEN_ROOT, "results.rds"))
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

ggplot2::ggsave(file.path("simulations", "all_plots", "6_segments_vs_jump.pdf"), plot,
                width = 16, height = 7)
cat("written", file.path(EVEN_ROOT, "segments_vs_jump.pdf"), "\n")
