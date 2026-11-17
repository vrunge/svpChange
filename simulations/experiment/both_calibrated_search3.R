## Dense (A, B) grid search around the ORIGINAL flat-threshold curve fit
## (A=1.683, B=6.209, C=4.939, from quantile_both_all_lengths.R), for the
## "both_calibrated" subtests mode
## (single-anchored FOCUS statistic, "both" pruning, threshold
## gamma(l) = A*log(l) + B*log(log(l)) + C evaluated at each candidate's
## OWN current length l). Many shapes, few replicates per shape ("much less
## replications, but much more simulations"): C is solved per (A,B) with a
## short, cheap bisection (not gridded -- it is almost purely a level shift
## once A, B and the 1% target are fixed), then the winner is re-validated
## with a high-precision bisection.
suppressPackageStartupMessages({
  library(svpChange)
  library(dplyr)
  library(ggplot2)
})

ROOT <- "/Users/vrunge/Library/CloudStorage/Dropbox/C_R_packages/svpChange"
setwd(ROOT)
source(file.path("simulations", "power_common.R"))

EXPERIMENT_ROOT <- file.path("simulations", "experiment")
N <- 1000L
TOLERANCE <- round(N * 0.0025)
FPR_TARGET <- 0.01
F1_JUMP <- 0.6
WORKERS <- power_default_workers()

FITTED_A <- 1.683
FITTED_B <- 6.209
FITTED_C <- 4.939

empirical_fpr <- function(A, B, C, reps, seed = 90001L) {
  hits <- power_lapply(seq_len(reps), function(i) {
    set_power_seed(seed + i)
    y <- rnorm(N)
    fit <- SVP(y, 0, subtests = "both_calibrated", abc_A = A, abc_B = B, abc_C = C)
    as.integer(length(fit$changepoints) > 1L)
  }, workers = WORKERS)
  mean(unlist(hits))
}

mean_updown_f1 <- function(A, B, C, jump, reps, seed = 70001L) {
  scores <- power_lapply(seq_len(reps), function(i) {
    set_power_seed(seed + i)
    mu <- generate_power_signal(N, "updown", jump)
    y <- mu + rnorm(N)
    fit <- SVP(y, 0, subtests = "both_calibrated", abc_A = A, abc_B = B, abc_C = C)
    truth <- normalise_boundaries(which(diff(mu) != 0), N)
    est <- normalise_boundaries(fit$changepoints, N)
    unname(paper_metrics(truth, est, TOLERANCE)["F1"])
  }, workers = WORKERS)
  mean(unlist(scores))
}

cheap_bisect_c <- function(A, B, C_lo = 8, C_hi = 22, target = FPR_TARGET,
                           reps = 80L, max_iter = 6L, seed = 90001L) {
  fpr_lo <- empirical_fpr(A, B, C_lo, reps, seed)
  fpr_hi <- empirical_fpr(A, B, C_hi, reps, seed)
  if (fpr_lo < target) return(list(C = C_lo, fpr = fpr_lo))
  if (fpr_hi > target) return(list(C = C_hi, fpr = fpr_hi))
  best <- list(C = C_hi, fpr = fpr_hi)
  for (iter in seq_len(max_iter)) {
    C_mid <- (C_lo + C_hi) / 2
    fpr_mid <- empirical_fpr(A, B, C_mid, reps, seed)
    if (abs(fpr_mid - target) < abs(best$fpr - target)) best <- list(C = C_mid, fpr = fpr_mid)
    if (fpr_mid > target) C_lo <- C_mid else C_hi <- C_mid
  }
  best
}

cat("=== Dense (A, B) grid around the fitted shape (A=1.683, B=6.209) ===\n")
A_grid <- seq(0.2, 3.2, by = 0.1)
B_grid <- seq(4.7, 7.7, by = 0.1)
grid <- expand.grid(A = A_grid, B = B_grid)
cat(sprintf("Grid size: %d (A,B) shapes\n", nrow(grid)))

t0 <- Sys.time()
results <- lapply(seq_len(nrow(grid)), function(i) {
  A <- grid$A[i]; B <- grid$B[i]
  cal <- cheap_bisect_c(A, B, reps = 80L)
  f1 <- mean_updown_f1(A, B, cal$C, F1_JUMP, reps = 80L)
  if (i %% 100 == 0) cat(sprintf("  ... %d / %d shapes done\n", i, nrow(grid)))
  c(A = A, B = B, C = cal$C, fpr = cal$fpr, f1 = f1)
})
cat(sprintf("Grid search took %.1f s\n", as.numeric(Sys.time() - t0, units = "secs")))

results <- as.data.frame(do.call(rbind, results))
saveRDS(results, file.path(EXPERIMENT_ROOT, "both_calibrated_search3_grid.rds"))

cat("\n=== Top 15 by F1, restricted to FPR in [0.3%, 2.5%] (noise band for 80 reps) ===\n")
plausible <- results[results$fpr >= 0.003 & results$fpr <= 0.025, ]
plausible <- plausible[order(-plausible$f1), ]
print(head(plausible, 15), row.names = FALSE)

winner <- plausible[1, ]
cat(sprintf("\nProvisional winner: A=%.2f B=%.2f C=%.3f (coarse FPR=%.3f, F1@0.6=%.4f)\n",
            winner$A, winner$B, winner$C, winner$fpr, winner$f1))

cat("\n=== High-precision re-calibration of the winning shape ===\n")
precise_bisect_c <- function(A, B, C_lo, C_hi, target = FPR_TARGET,
                             reps = 3000L, tol = 0.008, max_iter = 12L, seed = 90001L) {
  fpr_lo <- empirical_fpr(A, B, C_lo, reps, seed)
  fpr_hi <- empirical_fpr(A, B, C_hi, reps, seed)
  cat(sprintf("bracket: C=%.3f->FPR=%.4f | C=%.3f->FPR=%.4f\n", C_lo, fpr_lo, C_hi, fpr_hi))
  best <- list(C = C_hi, fpr = fpr_hi)
  for (iter in seq_len(max_iter)) {
    C_mid <- (C_lo + C_hi) / 2
    fpr_mid <- empirical_fpr(A, B, C_mid, reps, seed)
    cat(sprintf("iter %2d: C=%.4f -> FPR=%.4f\n", iter, C_mid, fpr_mid))
    if (abs(fpr_mid - target) < abs(best$fpr - target)) best <- list(C = C_mid, fpr = fpr_mid)
    if (abs(fpr_mid - target) < tol) break
    if (fpr_mid > target) C_lo <- C_mid else C_hi <- C_mid
  }
  best
}
precise <- precise_bisect_c(winner$A, winner$B, winner$C - 1.5, winner$C + 1.5)
final_f1 <- mean_updown_f1(winner$A, winner$B, precise$C, F1_JUMP, reps = 1000L)
best <- list(A = winner$A, B = winner$B, C = precise$C, fpr = precise$fpr, f1 = final_f1)
cat(sprintf("\nFinal winner: A=%.3f B=%.3f C=%.4f FPR=%.4f F1@0.6=%.4f (1000 reps)\n",
            best$A, best$B, best$C, best$fpr, best$f1))
saveRDS(best, file.path(EXPERIMENT_ROOT, "both_calibrated_best3.rds"))

cat("\n=== Final validation (5000 reps) ===\n")
final_fpr <- empirical_fpr(best$A, best$B, best$C, reps = 5000L, seed = 500001L)
cat(sprintf("Winner FPR at n=1000, 5000 reps: %.4f\n", final_fpr))
best$final_fpr_5000 <- final_fpr
saveRDS(best, file.path(EXPERIMENT_ROOT, "both_calibrated_best3.rds"))

cat("\n=== Full F1-vs-jump curve (updown) vs previous best2/best ===\n")
JUMP_SIZES <- seq(0.1, 2, 0.1)
prev_best <- tryCatch(readRDS(file.path(EXPERIMENT_ROOT, "both_calibrated_best.rds")),
                      error = function(e) NULL)

f1_curve <- function(A, B, C, label, reps = 300L) {
  vals <- vapply(JUMP_SIZES, function(j) mean_updown_f1(A, B, C, j, reps = reps),
                 numeric(1))
  data.frame(jump = JUMP_SIZES, F1 = vals, method = label)
}

curve_new <- f1_curve(best$A, best$B, best$C,
                      sprintf("dense grid near fitted shape (A=%.3f,B=%.3f,C=%.3f, FPR=%.1f%%)",
                              best$A, best$B, best$C, 100 * final_fpr))
curves <- curve_new
if (!is.null(prev_best)) {
  curve_prev <- f1_curve(prev_best$A, prev_best$B, prev_best$C,
                         sprintf("previous 4-shape search (A=%.3f,B=%.3f,C=%.3f, FPR~%.0f%%)",
                                 prev_best$A, prev_best$B, prev_best$C, 100 * prev_best$fpr))
  curves <- rbind(curve_prev, curve_new)
}
saveRDS(curves, file.path(EXPERIMENT_ROOT, "both_calibrated_search3_f1.rds"))

f1_plot <- ggplot2::ggplot(curves, ggplot2::aes(jump, F1, colour = method)) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::geom_point(size = 1.5) +
  ggplot2::labs(x = "Jump size", y = "F1 (updown scenario)",
                title = "both_calibrated: dense (A,B) grid near the fitted shape") +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "bottom", legend.direction = "vertical")
ggplot2::ggsave(file.path(EXPERIMENT_ROOT, "both_calibrated_search3_comparison.pdf"),
                f1_plot, width = 7, height = 6)

cat("\nDone. Outputs written to simulations/experiment/:\n")
cat("  both_calibrated_search3_grid.rds (full grid)\n")
cat("  both_calibrated_best3.rds\n")
cat("  both_calibrated_search3_f1.rds\n")
cat("  both_calibrated_search3_comparison.pdf\n")
