## Calibration for "both_mosum_sqrt": the MOSUM (moving-sum) validity test with
## a multiscale, length-adaptive penalty (Eichinger & Kirch 2018, Bernoulli
## 24(1), Theorem 2.1 asymptotic boundary, evaluated at the candidate's own
## length). gamma is a single additive correction on top of the exact
## theoretical boundary (gamma=0 reproduces it exactly); found by directly
## bisecting the empirical null rate to 1% under "both" pruning, exactly
## like every other method in this document.
source(file.path("simulations", "power_gaussian", "calibrate_modes.R"))
library(svpChange)
source(file.path("simulations", "power_common.R"))

N_CAL <- 1000L
TOLERANCE <- round(N_CAL * 0.0025)
ALPHA_TARGET <- 0.01

workers <- gaussian_calibration_workers()
if (workers < 2L) workers <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
cat("workers:", workers, "\n")

empirical_fpr <- function(gamma, reps, seed) {
  hits <- calibration_lapply(seq_len(reps), function(i) {
    set_power_seed(seed + i)
    y <- rnorm(N_CAL)
    fit <- SVP(y, gamma, subtests = "both_mosum_sqrt")
    as.integer(length(fit$changepoints) > 1L)
  }, workers)
  mean(unlist(hits))
}

mean_updown_f1 <- function(gamma, jump, reps, seed) {
  scores <- power_lapply(seq_len(reps), function(i) {
    set_power_seed(seed + i)
    mu <- generate_power_signal(N_CAL, "updown", jump)
    y <- mu + rnorm(N_CAL)
    fit <- SVP(y, gamma, subtests = "both_mosum_sqrt")
    truth <- normalise_boundaries(which(diff(mu) != 0), N_CAL)
    est <- normalise_boundaries(fit$changepoints, N_CAL)
    unname(paper_metrics(truth, est, TOLERANCE)["F1"])
  }, workers = workers)
  mean(unlist(scores))
}

cat("=== Bracket search around gamma=0 ===\n")
g_lo <- -1.0; g_hi <- 1.0
reps <- 3000L; seed <- 760000L
fpr_lo <- empirical_fpr(g_lo, reps, seed)
fpr_hi <- empirical_fpr(g_hi, reps, seed)
cat(sprintf("gamma=%.3f -> FPR=%.4f | gamma=%.3f -> FPR=%.4f\n", g_lo, fpr_lo, g_hi, fpr_hi))
while (fpr_lo < ALPHA_TARGET && g_lo > -10) {
  g_lo <- g_lo - 1.0
  fpr_lo <- empirical_fpr(g_lo, reps, seed)
  cat(sprintf("widening low: gamma=%.3f -> FPR=%.4f\n", g_lo, fpr_lo))
}
while (fpr_hi > ALPHA_TARGET && g_hi < 10) {
  g_hi <- g_hi + 1.0
  fpr_hi <- empirical_fpr(g_hi, reps, seed)
  cat(sprintf("widening high: gamma=%.3f -> FPR=%.4f\n", g_hi, fpr_hi))
}

cat("\n=== Bisecting gamma for a true 1% empirical FPR ===\n")
best <- list(gamma = g_hi, fpr = fpr_hi)
for (iter in 1:16) {
  g_mid <- (g_lo + g_hi) / 2
  fpr_mid <- empirical_fpr(g_mid, reps, seed)
  cat(sprintf("iter %2d: gamma=%.4f -> FPR=%.4f\n", iter, g_mid, fpr_mid))
  if (abs(fpr_mid - ALPHA_TARGET) < abs(best$fpr - ALPHA_TARGET)) best <- list(gamma = g_mid, fpr = fpr_mid)
  if (abs(fpr_mid - ALPHA_TARGET) < 0.0008) break
  if (fpr_mid > ALPHA_TARGET) g_lo <- g_mid else g_hi <- g_mid
}
gamma_star <- best$gamma
cat(sprintf("\nSelected gamma = %.4f (bisection FPR estimate = %.4f)\n", gamma_star, best$fpr))

cat("\n=== High-precision final validation (8000 reps) ===\n")
final_fpr <- empirical_fpr(gamma_star, 8000L, 970000L)
cat(sprintf("Final empirical FPR at n=1000 (8000 reps): %.4f\n", final_fpr))

cat("\n=== F1 at jump=0.6 on updown (2000 reps) ===\n")
final_f1 <- mean_updown_f1(gamma_star, jump = 0.6, reps = 2000L, seed = 1090000L)
cat(sprintf("F1@0.6 = %.4f\n", final_f1))

result <- list(gamma = gamma_star, final_fpr = final_fpr, final_f1_at_06 = final_f1)
saveRDS(result, file.path("simulations", "power_gaussian", "mosum_sqrt_calibration.rds"))
print(result)
