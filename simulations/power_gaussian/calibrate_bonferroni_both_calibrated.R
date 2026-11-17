## "both_calibrated" with a Bonferroni-corrected C, keeping the ORIGINAL
## fitted (A, B) shape fixed (A=1.683, B=6.209). The raw reuse
## of that fit (C=4.939) in both_calibrated gives ~75% empirical null FPR
## under "both" pruning, because gamma(l) = A*log(l)+B*log(log(l))+C is
## applied once per distinct candidate length l encountered in a run, not
## once per total-length evaluation. Treating this as K implicit per-length
## tests, a Bonferroni
## correction re-targets each one at alpha/K so the family-wise FPR over the
## whole run stays near alpha=1%: re-derive C (only C, A and B held fixed)
## from the (1 - alpha/K) quantile of the aggregate "both" statistic at
## n=1000, obtained via a Gumbel tail extrapolation (the quantile itself is
## too extreme to estimate directly from raw order statistics without an
## infeasible number of replicates).
source(file.path("simulations", "power_gaussian", "calibrate_modes.R"))
library(svpChange)

FITTED_A <- 1.683
FITTED_B <- 6.209
N_CAL <- 1000L
K_BONFERRONI <- 1000L          # one implicit test per candidate length l in [1, n]
ALPHA_TARGET <- 0.01
ALPHA_BONF <- ALPHA_TARGET / K_BONFERRONI

cat(sprintf("K = %d, alpha_target = %.4f, alpha_bonf = %.2e\n",
            K_BONFERRONI, ALPHA_TARGET, ALPHA_BONF))

compile_gaussian_critical_values()
workers <- gaussian_calibration_workers()
if (workers < 2L) {
  workers <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
}
cat("workers:", workers, "\n")

REPS <- 150000L
SEED <- 4210000L

cat("Simulating", REPS, "aggregate 'both' null maxima at n =", N_CAL, "...\n")
t0 <- Sys.time()
vals <- calibration_lapply(seq_len(REPS), function(i) {
  y <- simulate_gaussian_null(N_CAL, SEED + i)
  gaussian_critical_values_cpp(y)[["both"]]
}, workers)
vals <- unlist(vals, use.names = FALSE)
cat("done in", as.numeric(Sys.time() - t0, units = "secs"), "s\n")

saveRDS(vals, file.path("simulations", "power_gaussian",
                        "bonferroni_both_null_sample.rds"))

## --- Gumbel tail fit and extrapolation -----------------------------------

fit_gumbel_tail <- function(vals, tail_frac = 0.05) {
  N <- length(vals)
  sorted <- sort(vals, decreasing = TRUE) # sorted[1] = max
  k <- seq_len(floor(N * tail_frac))
  x <- sorted[k]
  S <- (k - 0.5) / N                     # empirical survival prob P(X > x)
  y <- -log(-log(1 - S))
  fit <- lm(x ~ y)
  list(fit = fit, a = unname(coef(fit)[2]), b = unname(coef(fit)[1]),
       r2 = summary(fit)$r.squared, tail_frac = tail_frac, n_tail = length(k))
}

gumbel <- fit_gumbel_tail(vals)
cat(sprintf("Gumbel tail fit: a = %.4f, b = %.4f, R^2 = %.6f (n_tail = %d)\n",
            gumbel$a, gumbel$b, gumbel$r2, gumbel$n_tail))

gumbel_quantile <- function(gumbel, alpha) {
  y_target <- -log(-log(1 - alpha))
  gumbel$b + gumbel$a * y_target
}

## Sanity check: compare the Gumbel-extrapolated 1% quantile to the DIRECT
## empirical 99% quantile from the same sample (should closely agree).
q_gumbel_01 <- gumbel_quantile(gumbel, 0.01)
q_empirical_01 <- as.numeric(quantile(vals, 0.99, type = 1))
cat(sprintf("Check at alpha=1%%: Gumbel-extrapolated = %.4f, direct empirical = %.4f (diff = %.4f)\n",
            q_gumbel_01, q_empirical_01, q_gumbel_01 - q_empirical_01))

gamma_bonf <- gumbel_quantile(gumbel, ALPHA_BONF)
cat(sprintf("Gumbel-extrapolated (1 - alpha_bonf) quantile at n=1000: %.4f\n", gamma_bonf))

C_new <- gamma_bonf - FITTED_A * log(N_CAL) - FITTED_B * log(log(N_CAL))
cat(sprintf("New C (A, B fixed at %.3f, %.3f): C = %.4f\n", FITTED_A, FITTED_B, C_new))

saveRDS(list(A = FITTED_A, B = FITTED_B, C = C_new, K = K_BONFERRONI,
            alpha_bonf = ALPHA_BONF, gamma_bonf = gamma_bonf, gumbel = gumbel),
        file.path("simulations", "power_gaussian", "bonferroni_both_calibrated.rds"))

## --- Empirical validation on the real discrete statistic -----------------

cat("\n=== Empirical validation of both_calibrated with Bonferroni C ===\n")
source(file.path("simulations", "power_common.R"))
fp <- 0L; reps_val <- 5000L
for (i in seq_len(reps_val)) {
  set_power_seed(900000L + i)
  y <- rnorm(N_CAL)
  fit <- SVP(y, 0, subtests = "both_calibrated",
            abc_A = FITTED_A, abc_B = FITTED_B, abc_C = C_new)
  if (length(fit$changepoints) > 1L) fp <- fp + 1L
}
cat(sprintf("Empirical FPR at n=1000 (%d reps): %.4f\n", reps_val, fp / reps_val))
