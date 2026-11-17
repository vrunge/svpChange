## Calibrate "both_rmin" (minimum-arm restricted CUSUM) to a true 1% null
## FPR at n=1000: search over r_min, bisecting gamma (direct Z^2-scale
## threshold) for each, then pick the r_min with the best updown F1 at
## jump=0.6 among those properly calibrated.
library(svpChange)
source(file.path("simulations", "power_common.R"))
source(file.path("simulations", "gaussian_common.R"))

N <- 1000L
TOLERANCE <- round(N * 0.0025)
ALPHA <- 0.01
WORKERS <- power_default_workers()

empirical_fpr <- function(r_min, gamma, reps, seed) {
  hits <- power_lapply(seq_len(reps), function(i) {
    set_power_seed(seed + i)
    y <- rnorm(N)
    fit <- SVP(y, gamma, subtests = "both_rmin", r_min = r_min)
    as.integer(length(fit$changepoints) > 1L)
  }, workers = WORKERS)
  mean(unlist(hits))
}

mean_updown_f1 <- function(r_min, gamma, jump, reps, seed) {
  scores <- power_lapply(seq_len(reps), function(i) {
    set_power_seed(seed + i)
    mu <- generate_power_signal(N, "updown", jump)
    y <- mu + rnorm(N)
    fit <- SVP(y, gamma, subtests = "both_rmin", r_min = r_min)
    truth <- normalise_boundaries(which(diff(mu) != 0), N)
    est <- normalise_boundaries(fit$changepoints, N)
    unname(paper_metrics(truth, est, TOLERANCE)["F1"])
  }, workers = WORKERS)
  mean(unlist(scores))
}

bisect_gamma <- function(r_min, g_lo, g_hi, target = ALPHA, reps = 1500L,
                         tol = 0.0015, max_iter = 16L, seed = 610000L, tag = "") {
  fpr_lo <- empirical_fpr(r_min, g_lo, reps, seed)
  fpr_hi <- empirical_fpr(r_min, g_hi, reps, seed)
  cat(sprintf("[%s] bracket: g=%.3f->FPR=%.4f | g=%.3f->FPR=%.4f\n",
              tag, g_lo, fpr_lo, g_hi, fpr_hi))
  best <- list(gamma = g_hi, fpr = fpr_hi)
  for (iter in seq_len(max_iter)) {
    g_mid <- (g_lo + g_hi) / 2
    fpr_mid <- empirical_fpr(r_min, g_mid, reps, seed)
    if (abs(fpr_mid - target) < abs(best$fpr - target)) best <- list(gamma = g_mid, fpr = fpr_mid)
    if (abs(fpr_mid - target) < tol) break
    if (fpr_mid > target) g_lo <- g_mid else g_hi <- g_mid
  }
  cat(sprintf("[%s] calibrated gamma=%.4f FPR=%.4f\n", tag, best$gamma, best$fpr))
  best
}

cat("=== Stage 1: bisect gamma for each r_min ===\n")
r_min_grid <- c(10, 20, 30, 40, 60, 80)
stage1 <- lapply(r_min_grid, function(rm) {
  tag <- sprintf("r_min=%d", rm)
  cal <- bisect_gamma(rm, g_lo = 15, g_hi = 30, tag = tag)
  f1 <- mean_updown_f1(rm, cal$gamma, jump = 0.6, reps = 800L, seed = 710000L)
  cat(sprintf("[%s] F1@0.6 = %.4f\n", tag, f1))
  c(r_min = rm, gamma = cal$gamma, fpr = cal$fpr, f1 = f1)
})
stage1 <- as.data.frame(do.call(rbind, stage1))
cat("\n=== Stage 1 results ===\n")
print(stage1[order(-stage1$f1), ], row.names = FALSE)

winner <- stage1[order(-stage1$f1), ][1, ]
cat(sprintf("\nProvisional winner: r_min=%d gamma=%.4f FPR=%.4f F1@0.6=%.4f\n",
            winner$r_min, winner$gamma, winner$fpr, winner$f1))

cat("\n=== High-precision recalibration of the winner ===\n")
final <- bisect_gamma(winner$r_min, g_lo = winner$gamma - 2, g_hi = winner$gamma + 2,
                      reps = 4000L, tol = 0.0008, max_iter = 16L, seed = 810000L,
                      tag = sprintf("final r_min=%d", winner$r_min))
final_fpr <- empirical_fpr(winner$r_min, final$gamma, reps = 8000L, seed = 910000L)
final_f1 <- mean_updown_f1(winner$r_min, final$gamma, jump = 0.6, reps = 2000L, seed = 1010000L)
cat(sprintf("\nFinal: r_min=%d gamma=%.4f FPR(8000 reps)=%.4f F1@0.6(2000 reps)=%.4f\n",
            winner$r_min, final$gamma, final_fpr, final_f1))

result <- list(r_min = winner$r_min, gamma = final$gamma, fpr = final_fpr, f1_at_06 = final_f1,
              stage1 = stage1)
saveRDS(result, file.path("simulations", "power_gaussian", "rmin_calibration.rds"))
utils::write.csv(stage1, file.path("simulations", "power_gaussian", "rmin_calibration_search.csv"),
                 row.names = FALSE)
cat("Saved.\n")
