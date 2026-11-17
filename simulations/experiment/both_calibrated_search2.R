## Wider, two-stage search for a well-calibrated (A, B, C) triple for the
## "both_calibrated" subtests mode (single-anchored FOCUS statistic, "both"
## pruning, threshold gamma(l) = A*log(l) + B*log(log(l)) + C evaluated at
## each candidate's OWN current length l). Supersedes
## both_calibrated_search.R's 4-shape search with a systematic grid plus a
## refinement stage around the best candidates.
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

empirical_fpr <- function(A, B, C, reps = 1000L, seed = 90001L) {
  hits <- power_lapply(seq_len(reps), function(i) {
    set_power_seed(seed + i)
    y <- rnorm(N)
    fit <- SVP(y, 0, subtests = "both_calibrated", abc_A = A, abc_B = B, abc_C = C)
    as.integer(length(fit$changepoints) > 1L)
  }, workers = WORKERS)
  mean(unlist(hits))
}

mean_updown_f1 <- function(A, B, C, jump, reps = 300L, seed = 70001L) {
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

bisect_c <- function(A, B, C_lo, C_hi, target = FPR_TARGET, reps = 1000L,
                     tol = 0.02, max_iter = 12L, seed = 90001L, tag = "") {
  fpr_lo <- empirical_fpr(A, B, C_lo, reps, seed)
  fpr_hi <- empirical_fpr(A, B, C_hi, reps, seed)
  cat(sprintf("[%s] bracket A=%.3f B=%.3f: C=%.3f->FPR=%.3f | C=%.3f->FPR=%.3f\n",
              tag, A, B, C_lo, fpr_lo, C_hi, fpr_hi))
  tries <- 0L
  while (fpr_hi > target && tries < 6L) {
    C_hi <- C_hi + (C_hi - C_lo)
    fpr_hi <- empirical_fpr(A, B, C_hi, reps, seed)
    tries <- tries + 1L
    cat(sprintf("[%s] widened: C=%.3f->FPR=%.3f\n", tag, C_hi, fpr_hi))
  }
  if (fpr_lo < target) {
    return(list(A = A, B = B, C = C_lo, fpr = fpr_lo))
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

cat("=== Stage 1: coarse grid over (A, B) ===\n")
A_grid <- c(0, 1, 2, 3, 4)
B_grid <- c(-2, 0, 2, 4)
grid <- expand.grid(A = A_grid, B = B_grid)

stage1 <- lapply(seq_len(nrow(grid)), function(i) {
  A <- grid$A[i]; B <- grid$B[i]
  tag <- sprintf("A=%g,B=%g", A, B)
  cal <- bisect_c(A, B, C_lo = 2, C_hi = 30, reps = 1000L, tol = 0.02, tag = tag)
  f1 <- mean_updown_f1(cal$A, cal$B, cal$C, F1_JUMP, reps = 300L)
  cat(sprintf("[%s] calibrated C=%.3f FPR=%.3f F1@0.6=%.4f\n",
              tag, cal$C, cal$fpr, f1))
  c(A = cal$A, B = cal$B, C = cal$C, fpr = cal$fpr, f1 = f1)
})
stage1 <- as.data.frame(do.call(rbind, stage1))
saveRDS(stage1, file.path(EXPERIMENT_ROOT, "both_calibrated_search2_stage1.rds"))

cat("\n=== Stage 1 results (sorted by F1) ===\n")
print(stage1[order(-stage1$f1), ], row.names = FALSE)

cat("\n=== Stage 2: refine around the top 5 shapes ===\n")
top5 <- stage1[order(-stage1$f1), ][seq_len(min(5L, nrow(stage1))), ]

refine_offsets <- expand.grid(dA = c(-0.5, 0, 0.5), dB = c(-1, 0, 1))
stage2_candidates <- do.call(rbind, lapply(seq_len(nrow(top5)), function(i) {
  A0 <- top5$A[i]; B0 <- top5$B[i]
  data.frame(A = pmax(A0 + refine_offsets$dA, 0), B = B0 + refine_offsets$dB)
}))
stage2_candidates <- unique(stage2_candidates)

stage2 <- lapply(seq_len(nrow(stage2_candidates)), function(i) {
  A <- stage2_candidates$A[i]; B <- stage2_candidates$B[i]
  tag <- sprintf("refine A=%.2f,B=%.2f", A, B)
  cal <- bisect_c(A, B, C_lo = 2, C_hi = 30, reps = 3000L, tol = 0.008, tag = tag)
  f1 <- mean_updown_f1(cal$A, cal$B, cal$C, F1_JUMP, reps = 500L)
  cat(sprintf("[%s] calibrated C=%.3f FPR=%.3f F1@0.6=%.4f\n",
              tag, cal$C, cal$fpr, f1))
  c(A = cal$A, B = cal$B, C = cal$C, fpr = cal$fpr, f1 = f1)
})
stage2 <- as.data.frame(do.call(rbind, stage2))
saveRDS(stage2, file.path(EXPERIMENT_ROOT, "both_calibrated_search2_stage2.rds"))

cat("\n=== Stage 2 results (sorted by F1) ===\n")
print(stage2[order(-stage2$f1), ], row.names = FALSE)

best <- as.list(stage2[order(-stage2$f1), ][1, ])
cat(sprintf("\nBest overall: A=%.3f B=%.3f C=%.3f FPR=%.3f F1@0.6=%.4f\n",
            best$A, best$B, best$C, best$fpr, best$f1))
saveRDS(best, file.path(EXPERIMENT_ROOT, "both_calibrated_best2.rds"))

cat("\n=== Final validation of the winner (higher precision) ===\n")
final_fpr <- empirical_fpr(best$A, best$B, best$C, reps = 5000L, seed = 500001L)
cat(sprintf("Winner FPR at n=1000, 5000 reps: %.4f\n", final_fpr))
best$final_fpr_5000 <- final_fpr
saveRDS(best, file.path(EXPERIMENT_ROOT, "both_calibrated_best2.rds"))

cat("\n=== Full F1-vs-jump curve (updown) for the winner vs stage-1 previous best ===\n")
JUMP_SIZES <- seq(0.1, 2, 0.1)
prev_best <- readRDS(file.path(EXPERIMENT_ROOT, "both_calibrated_best.rds"))

f1_curve <- function(A, B, C, label) {
  vals <- vapply(JUMP_SIZES, function(j) mean_updown_f1(A, B, C, j, reps = 300L),
                 numeric(1))
  data.frame(jump = JUMP_SIZES, F1 = vals, method = label)
}

curve_prev <- f1_curve(prev_best$A, prev_best$B, prev_best$C,
                       sprintf("previous search (A=%.3f,B=%.3f,C=%.3f, FPR~%.0f%%)",
                               prev_best$A, prev_best$B, prev_best$C,
                               100 * prev_best$fpr))
curve_new <- f1_curve(best$A, best$B, best$C,
                      sprintf("wider search (A=%.3f,B=%.3f,C=%.3f, FPR=%.1f%%)",
                              best$A, best$B, best$C, 100 * final_fpr))

f1_data <- rbind(curve_prev, curve_new)
saveRDS(f1_data, file.path(EXPERIMENT_ROOT, "both_calibrated_search2_f1.rds"))

f1_plot <- ggplot2::ggplot(f1_data, ggplot2::aes(jump, F1, colour = method)) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::geom_point(size = 1.5) +
  ggplot2::labs(x = "Jump size", y = "F1 (updown scenario)",
                title = "both_calibrated: previous vs wider (A,B,C) search") +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "bottom", legend.direction = "vertical")
ggplot2::ggsave(file.path(EXPERIMENT_ROOT, "both_calibrated_search2_comparison.pdf"),
                f1_plot, width = 7, height = 6)

cat("\nDone. Outputs written to simulations/experiment/:\n")
cat("  both_calibrated_search2_stage1.rds\n")
cat("  both_calibrated_search2_stage2.rds\n")
cat("  both_calibrated_best2.rds\n")
cat("  both_calibrated_search2_f1.rds\n")
cat("  both_calibrated_search2_comparison.pdf\n")
