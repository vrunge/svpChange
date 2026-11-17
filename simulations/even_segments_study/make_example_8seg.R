## Example figure for the paper section: eight segments, small jump size
## 0.6. PELT infers 5 changepoints in total (3 correctly matching a true
## change, 2 spurious). SVP infers only 3 changepoints (2 correct, 1
## spurious) -- fewer detections overall than PELT. Regenerates the
## selected series from its simulation seed and plots the true change
## points, the PELT detections and the SVP detections.
##
## Run from the package root:
##   Rscript simulations/even_segments_study/make_example_8seg.R

source(file.path("simulations", "power_common.R"))
source(file.path("simulations", "even_segments_study", "run_study.R"))

EXAMPLE_SEED <- 7225162L

results <- readRDS(file.path(EVEN_ROOT, "results.rds"))
row <- results[results$seed == EXAMPLE_SEED, , drop = FALSE]
stopifnot(nrow(row) == 1L)
segments <- row$segments[[1L]]
jump <- row$jump[[1L]]
mu <- even_updown_signal(EVEN_N, segments, jump)
set_power_seed(EXAMPLE_SEED)
y <- mu + stats::rnorm(EVEN_N)
truth <- even_changepoints(EVEN_N, segments)
pelt_cps <- row$pelt_changepoints[[1L]]
svp_cps <- row$svp_changepoints[[1L]]

grDevices::pdf(file.path("simulations", "all_plots", "6_example_missed_change.pdf"),
               width = 11, height = 4.6)
graphics::par(mar = c(4, 4, 2.2, 1))
graphics::plot(
  seq_len(EVEN_N), y, pch = 19, cex = 0.35, col = "grey55",
  xlab = "Time", ylab = "Value",
  main = sprintf("Eight segments, jump size %.1f", jump)
)
graphics::lines(seq_len(EVEN_N), mu, col = "black", lwd = 1.6)
graphics::abline(v = truth, col = "grey25", lty = 3, lwd = 1.4)
graphics::abline(v = pelt_cps, col = "dodgerblue3", lty = 2, lwd = 4)
graphics::abline(v = svp_cps, col = "firebrick3", lty = 1, lwd = 1.5)
graphics::legend(
  "topright",
  legend = c("data", "true mean", "true changes", "PELT changes",
             "SVP changes"),
  col = c("grey55", "black", "grey25", "dodgerblue3", "firebrick3"),
  pch = c(19, NA, NA, NA, NA), lty = c(NA, 1, 3, 2, 1),
  lwd = c(NA, 1.6, 1.4, 4, 1.5), cex = 0.85, bg = "white"
)
grDevices::dev.off()

cat("truth:", truth, "\n")
cat("PELT :", pelt_cps, " F1", round(row$pelt_f1, 3), "\n")
cat("SVP  :", svp_cps, " F1", round(row$svp_f1, 3), "\n")
