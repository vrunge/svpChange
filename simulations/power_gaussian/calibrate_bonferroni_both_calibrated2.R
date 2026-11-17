## Refinement of calibrate_bonferroni_both_calibrated.R: instead of a fixed
## K=1000 Bonferroni divisor (alpha/1000), which overshot badly (0.16% FPR
## instead of 1%), solve directly for the effective per-length alpha that
## brings both_calibrated (A=1.683, B=6.209 fixed) back to a true 1% null
## FPR, via bisection on the real discrete statistic. The implied K is
## reported for reference using the Gumbel tail model already fit.
source(file.path("simulations", "power_gaussian", "calibrate_modes.R"))
library(svpChange)
source(file.path("simulations", "power_common.R"))

FITTED_A <- 1.683
FITTED_B <- 6.209
N_CAL <- 1000L
ALPHA_TARGET <- 0.01

gumbel_info <- readRDS(file.path("simulations", "power_gaussian",
                                 "bonferroni_both_calibrated.rds"))
gumbel <- gumbel_info$gumbel
cat(sprintf("Reusing Gumbel tail fit: a=%.4f, b=%.4f, R^2=%.6f\n",
            gumbel$a, gumbel$b, gumbel$r2))

workers <- gaussian_calibration_workers()
if (workers < 2L) workers <- max(1L, parallel::detectCores(logical = FALSE) - 1L)

empirical_fpr <- function(A, B, C, reps, seed) {
  hits <- calibration_lapply(seq_len(reps), function(i) {
    set_power_seed(seed + i)
    y <- rnorm(N_CAL)
    fit <- SVP(y, 0, subtests = "both_calibrated", abc_A = A, abc_B = B, abc_C = C)
    as.integer(length(fit$changepoints) > 1L)
  }, workers)
  mean(unlist(hits))
}

cat("=== Bisecting C for a true 1% empirical FPR (A, B fixed) ===\n")
C_lo <- 14.0; C_hi <- 16.0   # brackets around the earlier rough C=14.939 (~1.2% FPR)
reps <- 4000L; seed <- 950000L
fpr_lo <- empirical_fpr(FITTED_A, FITTED_B, C_lo, reps, seed)
fpr_hi <- empirical_fpr(FITTED_A, FITTED_B, C_hi, reps, seed)
cat(sprintf("bracket: C=%.3f -> FPR=%.4f | C=%.3f -> FPR=%.4f\n", C_lo, fpr_lo, C_hi, fpr_hi))
best <- list(C = C_hi, fpr = fpr_hi)
for (iter in 1:14) {
  C_mid <- (C_lo + C_hi) / 2
  fpr_mid <- empirical_fpr(FITTED_A, FITTED_B, C_mid, reps, seed)
  cat(sprintf("iter %2d: C=%.4f -> FPR=%.4f\n", iter, C_mid, fpr_mid))
  if (abs(fpr_mid - ALPHA_TARGET) < abs(best$fpr - ALPHA_TARGET)) best <- list(C = C_mid, fpr = fpr_mid)
  if (abs(fpr_mid - ALPHA_TARGET) < 0.0008) break
  if (fpr_mid > ALPHA_TARGET) C_lo <- C_mid else C_hi <- C_mid
}
C_star <- best$C
cat(sprintf("\nSelected C = %.4f (bisection FPR estimate = %.4f)\n", C_star, best$fpr))

cat("\n=== High-precision final validation (8000 reps) ===\n")
final_fpr <- empirical_fpr(FITTED_A, FITTED_B, C_star, 8000L, 960000L)
cat(sprintf("Final empirical FPR at n=1000 (8000 reps): %.4f\n", final_fpr))

## Back out the implied per-length alpha and K via the Gumbel tail model.
gamma_star <- C_star + FITTED_A * log(N_CAL) + FITTED_B * log(log(N_CAL))
y_val <- (gamma_star - gumbel$b) / gumbel$a
alpha_implied <- 1 - exp(-exp(-y_val))
K_implied <- ALPHA_TARGET / alpha_implied
cat(sprintf("\nImplied per-length alpha' = %.3e  (K' = alpha_target/alpha' = %.1f)\n",
            alpha_implied, K_implied))
cat("(K=1000 was the naive guess; the data say the effectively independent\n")
cat(" number of per-length tests is much smaller, consistent with strong\n")
cat(" correlation between tests at nearby candidate lengths.)\n")

result <- list(A = FITTED_A, B = FITTED_B, C = C_star,
              alpha_implied = alpha_implied, K_implied = K_implied,
              final_fpr = final_fpr, gamma_star = gamma_star)
saveRDS(result, file.path("simulations", "power_gaussian",
                          "bonferroni_both_calibrated_tight.rds"))
print(result)
