## Multiplicative-scale Bonferroni correction for "both_calibrated": instead
## of re-deriving the additive constant C, scale the WHOLE fitted threshold
## by a single factor w,
##
##   gamma(l) = w * (A*log(l) + B*log(log(l)) + C),   l = candidate's own
##   length, with the ORIGINAL fitted (A, B, C) = (1.683, 6.209, 4.939) held
##   fixed and unscaled individually.
##
## w is found by bisecting the true empirical null FPR of the real discrete
## statistic under "both" pruning at n=1000, directly (no Gumbel tail model,
## no independence/K assumption) -- the same direct-bisection methodology
## used throughout this document (both_rmin, the additive-C variant).
source(file.path("simulations", "power_gaussian", "calibrate_modes.R"))
library(svpChange)
source(file.path("simulations", "power_common.R"))
source(file.path("simulations", "gaussian_common.R"))

FITTED_A <- 1.683
FITTED_B <- 6.209
FITTED_C <- 4.939
N_CAL <- 1000L
TOLERANCE <- round(N_CAL * 0.0025)
ALPHA_TARGET <- 0.01

workers <- gaussian_calibration_workers()
if (workers < 2L) workers <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
cat("workers:", workers, "\n")

empirical_fpr <- function(w, reps, seed) {
  hits <- calibration_lapply(seq_len(reps), function(i) {
    set_power_seed(seed + i)
    y <- rnorm(N_CAL)
    fit <- SVP(y, 0, subtests = "both_calibrated",
              abc_A = w * FITTED_A, abc_B = w * FITTED_B, abc_C = w * FITTED_C)
    as.integer(length(fit$changepoints) > 1L)
  }, workers)
  mean(unlist(hits))
}

mean_updown_f1 <- function(w, jump, reps, seed) {
  scores <- power_lapply(seq_len(reps), function(i) {
    set_power_seed(seed + i)
    mu <- generate_power_signal(N_CAL, "updown", jump)
    y <- mu + rnorm(N_CAL)
    fit <- SVP(y, 0, subtests = "both_calibrated",
              abc_A = w * FITTED_A, abc_B = w * FITTED_B, abc_C = w * FITTED_C)
    truth <- normalise_boundaries(which(diff(mu) != 0), N_CAL)
    est <- normalise_boundaries(fit$changepoints, N_CAL)
    unname(paper_metrics(truth, est, TOLERANCE)["F1"])
  }, workers = workers)
  mean(unlist(scores))
}

cat("=== Bracket search ===\n")
w_lo <- 1.0; w_hi <- 2.0
reps <- 2000L; seed <- 720000L
fpr_lo <- empirical_fpr(w_lo, reps, seed)
fpr_hi <- empirical_fpr(w_hi, reps, seed)
cat(sprintf("w=%.3f -> FPR=%.4f | w=%.3f -> FPR=%.4f\n", w_lo, fpr_lo, w_hi, fpr_hi))
while (fpr_hi > ALPHA_TARGET && w_hi < 8) {
  w_hi <- w_hi * 1.5
  fpr_hi <- empirical_fpr(w_hi, reps, seed)
  cat(sprintf("widening: w=%.3f -> FPR=%.4f\n", w_hi, fpr_hi))
}

cat("\n=== Bisecting w for a true 1% empirical FPR ===\n")
best <- list(w = w_hi, fpr = fpr_hi)
for (iter in 1:16) {
  w_mid <- (w_lo + w_hi) / 2
  fpr_mid <- empirical_fpr(w_mid, reps, seed)
  cat(sprintf("iter %2d: w=%.4f -> FPR=%.4f\n", iter, w_mid, fpr_mid))
  if (abs(fpr_mid - ALPHA_TARGET) < abs(best$fpr - ALPHA_TARGET)) best <- list(w = w_mid, fpr = fpr_mid)
  if (abs(fpr_mid - ALPHA_TARGET) < 0.0008) break
  if (fpr_mid > ALPHA_TARGET) w_lo <- w_mid else w_hi <- w_mid
}
w_star <- best$w
cat(sprintf("\nSelected w = %.4f (bisection FPR estimate = %.4f)\n", w_star, best$fpr))

cat("\n=== High-precision final validation (8000 reps) ===\n")
final_fpr <- empirical_fpr(w_star, 8000L, 960000L)
cat(sprintf("Final empirical FPR at n=1000 (8000 reps): %.4f\n", final_fpr))

cat("\n=== F1 at jump=0.6 on updown (2000 reps) ===\n")
final_f1 <- mean_updown_f1(w_star, jump = 0.6, reps = 2000L, seed = 1010000L)
cat(sprintf("F1@0.6 = %.4f\n", final_f1))

result <- list(A = FITTED_A, B = FITTED_B, C = FITTED_C, w = w_star,
              final_fpr = final_fpr, final_f1_at_06 = final_f1)
saveRDS(result, file.path("simulations", "power_gaussian",
                          "bonferroni_both_calibrated_w.rds"))
print(result)
