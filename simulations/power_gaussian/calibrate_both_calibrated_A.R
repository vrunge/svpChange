## Recalibration of "both_calibrated" by tuning A DIRECTLY, instead of a
## global multiplicative scale w (calibrate_bonferroni_both_calibrated_w.R)
## or an additive-only C (calibrate_bonferroni_both_calibrated2.R). B and C
## are held at their original fitted values; only the log(l) coefficient A
## is re-derived,
##
##   gamma(l) = A * log(l) + B * log(log(l)) + C,   l = candidate's own
##   length, (B, C) = (6.209, 4.939) fixed at the original fit.
##
## A is found by bisecting the true empirical null FPR of the real discrete
## statistic under "both" pruning at n=1000, directly (same methodology as
## every other calibration in this document).
source(file.path("simulations", "power_gaussian", "calibrate_modes.R"))
library(svpChange)
source(file.path("simulations", "power_common.R"))
source(file.path("simulations", "gaussian_common.R"))

FITTED_B <- 6.209
FITTED_C <- 4.939
N_CAL <- 1000L
TOLERANCE <- round(N_CAL * 0.0025)
ALPHA_TARGET <- 0.01

workers <- gaussian_calibration_workers()
if (workers < 2L) workers <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
cat("workers:", workers, "\n")

empirical_fpr <- function(A, reps, seed) {
  hits <- calibration_lapply(seq_len(reps), function(i) {
    set_power_seed(seed + i)
    y <- rnorm(N_CAL)
    fit <- SVP(y, 0, subtests = "both_calibrated",
              abc_A = A, abc_B = FITTED_B, abc_C = FITTED_C)
    as.integer(length(fit$changepoints) > 1L)
  }, workers)
  mean(unlist(hits))
}

mean_updown_f1 <- function(A, jump, reps, seed) {
  scores <- power_lapply(seq_len(reps), function(i) {
    set_power_seed(seed + i)
    mu <- generate_power_signal(N_CAL, "updown", jump)
    y <- mu + rnorm(N_CAL)
    fit <- SVP(y, 0, subtests = "both_calibrated",
              abc_A = A, abc_B = FITTED_B, abc_C = FITTED_C)
    truth <- normalise_boundaries(which(diff(mu) != 0), N_CAL)
    est <- normalise_boundaries(fit$changepoints, N_CAL)
    unname(paper_metrics(truth, est, TOLERANCE)["F1"])
  }, workers = workers)
  mean(unlist(scores))
}

cat("=== Bracket search ===\n")
A_lo <- 1.683; A_hi <- 3.5
reps <- 2000L; seed <- 730000L
fpr_lo <- empirical_fpr(A_lo, reps, seed)
fpr_hi <- empirical_fpr(A_hi, reps, seed)
cat(sprintf("A=%.3f -> FPR=%.4f | A=%.3f -> FPR=%.4f\n", A_lo, fpr_lo, A_hi, fpr_hi))
while (fpr_hi > ALPHA_TARGET && A_hi < 20) {
  A_hi <- A_hi * 1.5
  fpr_hi <- empirical_fpr(A_hi, reps, seed)
  cat(sprintf("widening: A=%.3f -> FPR=%.4f\n", A_hi, fpr_hi))
}

cat("\n=== Bisecting A for a true 1% empirical FPR ===\n")
best <- list(A = A_hi, fpr = fpr_hi)
for (iter in 1:16) {
  A_mid <- (A_lo + A_hi) / 2
  fpr_mid <- empirical_fpr(A_mid, reps, seed)
  cat(sprintf("iter %2d: A=%.4f -> FPR=%.4f\n", iter, A_mid, fpr_mid))
  if (abs(fpr_mid - ALPHA_TARGET) < abs(best$fpr - ALPHA_TARGET)) best <- list(A = A_mid, fpr = fpr_mid)
  if (abs(fpr_mid - ALPHA_TARGET) < 0.0008) break
  if (fpr_mid > ALPHA_TARGET) A_lo <- A_mid else A_hi <- A_mid
}
A_star <- best$A
cat(sprintf("\nSelected A = %.4f (bisection FPR estimate = %.4f)\n", A_star, best$fpr))

cat("\n=== High-precision final validation (8000 reps) ===\n")
final_fpr <- empirical_fpr(A_star, 8000L, 960000L)
cat(sprintf("Final empirical FPR at n=1000 (8000 reps): %.4f\n", final_fpr))

cat("\n=== F1 at jump=0.6 on updown (2000 reps) ===\n")
final_f1 <- mean_updown_f1(A_star, jump = 0.6, reps = 2000L, seed = 1010000L)
cat(sprintf("F1@0.6 = %.4f\n", final_f1))

result <- list(A = A_star, B = FITTED_B, C = FITTED_C,
              final_fpr = final_fpr, final_f1_at_06 = final_f1)
saveRDS(result, file.path("simulations", "power_gaussian",
                          "both_calibrated_A_calibration.rds"))
print(result)
