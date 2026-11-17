## Search for a well-calibrated (A, B, C) triple for the "both_calibrated"
## subtests mode (single-anchored FOCUS statistic, "both" pruning, threshold
## gamma(l) = A*log(l) + B*log(log(l)) + C evaluated at each candidate's OWN
## current length l). Produces two curves:
##   1. "fitted"    A=1.683, B=6.209, C=4.939 (the flat-threshold curve fit
##                  from quantile_both_all_lengths.R, known to be badly
##                  miscalibrated when reused at each candidate's own length
##                  under "both" pruning: ~76.7%
##                  empirical FPR at n=1000 instead of 1%).
##   2. "recalibrated" a triple found by grid-searching (A,B) shapes and
##                  bisecting C on each shape to hit ~1% empirical FPR at
##                  n=1000, then picking the shape with the best mean F1 on
##                  the updown scenario.
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

empirical_fpr <- function(A, B, C, reps = 500L, seed = 90001L) {
  hits <- power_lapply(seq_len(reps), function(i) {
    set_power_seed(seed + i)
    y <- rnorm(N)
    fit <- SVP(y, 0, subtests = "both_calibrated", abc_A = A, abc_B = B, abc_C = C)
    as.integer(length(fit$changepoints) > 1L)
  }, workers = WORKERS)
  mean(unlist(hits))
}

mean_updown_f1 <- function(A, B, C, jump, reps = 150L, seed = 70001L) {
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

bisect_c <- function(A, B, C_lo, C_hi, target = FPR_TARGET, reps = 500L,
                     tol = 0.03, max_iter = 12L, seed = 90001L, tag = "") {
  fpr_lo <- empirical_fpr(A, B, C_lo, reps, seed)
  fpr_hi <- empirical_fpr(A, B, C_hi, reps, seed)
  cat(sprintf("[%s] bracket A=%.3f B=%.3f: C=%.3f->FPR=%.3f | C=%.3f->FPR=%.3f\n",
              tag, A, B, C_lo, fpr_lo, C_hi, fpr_hi))
  if (fpr_lo < target) {
    cat(sprintf("[%s] lower bound already under target; using C_lo\n", tag))
    return(list(A = A, B = B, C = C_lo, fpr = fpr_lo))
  }
  if (fpr_hi > target) {
    cat(sprintf("[%s] upper bound still over target; widening\n", tag))
    C_hi <- C_hi + (C_hi - C_lo)
    fpr_hi <- empirical_fpr(A, B, C_hi, reps, seed)
    cat(sprintf("[%s] widened: C=%.3f->FPR=%.3f\n", tag, C_hi, fpr_hi))
  }
  best <- list(A = A, B = B, C = C_hi, fpr = fpr_hi)
  for (iter in seq_len(max_iter)) {
    C_mid <- (C_lo + C_hi) / 2
    fpr_mid <- empirical_fpr(A, B, C_mid, reps, seed)
    cat(sprintf("[%s] iter %2d: C=%.4f -> FPR=%.3f\n", tag, iter, C_mid, fpr_mid))
    if (abs(fpr_mid - target) < abs(best$fpr - target)) {
      best <- list(A = A, B = B, C = C_mid, fpr = fpr_mid)
    }
    if (abs(fpr_mid - target) < tol) break
    if (fpr_mid > target) C_lo <- C_mid else C_hi <- C_mid
  }
  best
}

cat("=== Calibrating candidate shapes (A, B) via bisection on C ===\n")
shapes <- list(
  fitted_shape = c(A = FITTED_A, B = FITTED_B),
  flat         = c(A = 0,        B = 0),
  mild_growth  = c(A = 1,        B = 1),
  log_only     = c(A = 2,        B = 0)
)

calibrated <- lapply(names(shapes), function(nm) {
  A <- shapes[[nm]][["A"]]; B <- shapes[[nm]][["B"]]
  bisect_c(A, B, C_lo = FITTED_C, C_hi = FITTED_C + 20, tag = nm)
})
names(calibrated) <- names(shapes)

cat("\n=== Calibrated shapes ===\n")
for (nm in names(calibrated)) {
  cc <- calibrated[[nm]]
  cat(sprintf("%-14s A=%.3f B=%.3f C=%.3f -> FPR=%.3f\n",
              nm, cc$A, cc$B, cc$C, cc$fpr))
}

cat(sprintf("\n=== Scoring shapes on updown F1 at jump=%.1f ===\n", F1_JUMP))
for (nm in names(calibrated)) {
  cc <- calibrated[[nm]]
  f1 <- mean_updown_f1(cc$A, cc$B, cc$C, F1_JUMP)
  calibrated[[nm]]$f1 <- f1
  cat(sprintf("%-14s F1=%.4f\n", nm, f1))
}

best_name <- names(calibrated)[which.max(vapply(calibrated, `[[`, numeric(1), "f1"))]
best <- calibrated[[best_name]]
cat(sprintf("\nBest shape: %s (A=%.3f, B=%.3f, C=%.3f, FPR=%.3f, F1@0.6=%.4f)\n",
            best_name, best$A, best$B, best$C, best$fpr, best$f1))

saveRDS(calibrated, file.path(EXPERIMENT_ROOT, "both_calibrated_search_shapes.rds"))
saveRDS(best, file.path(EXPERIMENT_ROOT, "both_calibrated_best.rds"))

cat("\n=== Building full F1-vs-jump curves (updown) for both curves ===\n")
JUMP_SIZES <- seq(0.1, 2, 0.1)

f1_curve <- function(A, B, C, label) {
  vals <- vapply(JUMP_SIZES, function(j) mean_updown_f1(A, B, C, j, reps = 150L),
                 numeric(1))
  data.frame(jump = JUMP_SIZES, F1 = vals, method = label)
}

curve_fitted <- f1_curve(FITTED_A, FITTED_B, FITTED_C,
                         sprintf("SVP both_calibrated fitted (A=%.3f,B=%.3f,C=%.3f, FPR=%.0f%%)",
                                 FITTED_A, FITTED_B, FITTED_C,
                                 100 * calibrated$fitted_shape$fpr))
curve_best <- f1_curve(best$A, best$B, best$C,
                       sprintf("SVP both_calibrated recalibrated (A=%.3f,B=%.3f,C=%.3f, FPR=%.0f%%)",
                               best$A, best$B, best$C, 100 * best$fpr))

f1_data <- rbind(curve_fitted, curve_best)
saveRDS(f1_data, file.path(EXPERIMENT_ROOT, "both_calibrated_f1_comparison.rds"))

f1_plot <- ggplot2::ggplot(f1_data, ggplot2::aes(jump, F1, colour = method)) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::geom_point(size = 1.5) +
  ggplot2::labs(x = "Jump size", y = "F1 (updown scenario)",
                title = "both_calibrated: fitted vs recalibrated (A,B,C)") +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "bottom", legend.direction = "vertical")
ggplot2::ggsave(file.path(EXPERIMENT_ROOT, "both_calibrated_ABC_comparison.pdf"),
                f1_plot, width = 7, height = 6)

cat("\n=== Building gamma(l) vs l curves for both parameter sets ===\n")
gamma_of_l <- function(A, B, C, l) A * log(l) + B * log(log(l)) + C
l_grid <- seq(10, 1000, 10)
gamma_data <- rbind(
  data.frame(l = l_grid,
             gamma = gamma_of_l(FITTED_A, FITTED_B, FITTED_C, l_grid),
             method = "fitted (A=1.683,B=6.209,C=4.939)"),
  data.frame(l = l_grid,
             gamma = gamma_of_l(best$A, best$B, best$C, l_grid),
             method = sprintf("recalibrated (A=%.3f,B=%.3f,C=%.3f)",
                              best$A, best$B, best$C))
)
saveRDS(gamma_data, file.path(EXPERIMENT_ROOT, "both_calibrated_gamma_comparison.rds"))

gamma_plot <- ggplot2::ggplot(gamma_data, ggplot2::aes(l, gamma, colour = method)) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::labs(x = "Candidate length l", y = "gamma(l)",
                title = "both_calibrated threshold gamma(l): fitted vs recalibrated") +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "bottom", legend.direction = "vertical")
ggplot2::ggsave(file.path(EXPERIMENT_ROOT, "both_calibrated_gamma_ABC_comparison.pdf"),
                gamma_plot, width = 7, height = 6)

cat("\nDone. Outputs written to simulations/experiment/:\n")
cat("  both_calibrated_search_shapes.rds\n")
cat("  both_calibrated_best.rds\n")
cat("  both_calibrated_f1_comparison.rds\n")
cat("  both_calibrated_ABC_comparison.pdf\n")
cat("  both_calibrated_gamma_comparison.rds\n")
cat("  both_calibrated_gamma_ABC_comparison.pdf\n")
