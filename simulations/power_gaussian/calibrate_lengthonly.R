## Calibration for "both_lengthonly" (SVP_validity_rules.tex proposal):
## f_Gauss(y_{s:t}) <= [q + sqrt(2*log(e*n/m))]^2, m = candidate's own
## length, n = total series length, dispatched with "both" pruning (a
## heuristic shortcut to this raw statistic, like every other both-family
## method here). q found by bisecting the true empirical FPR of the real
## discrete statistic directly against subtests = "both_lengthonly".
source(file.path("simulations", "power_gaussian", "calibrate_modes.R"))
library(svpChange)
source(file.path("simulations", "power_common.R"))

N_CAL <- 1000L
TOLERANCE <- round(N_CAL * 0.0025)
ALPHA_TARGET <- 0.01

workers <- gaussian_calibration_workers()
if (workers < 2L) workers <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
cat("workers:", workers, "\n")

empirical_fpr <- function(q, reps, seed) {
  hits <- calibration_lapply(seq_len(reps), function(i) {
    set_power_seed(seed + i)
    y <- rnorm(N_CAL)
    fit <- SVP(y, q, subtests = "both_lengthonly")
    as.integer(length(fit$changepoints) > 1L)
  }, workers)
  mean(unlist(hits))
}

mean_updown_f1 <- function(q, jump, reps, seed) {
  scores <- power_lapply(seq_len(reps), function(i) {
    set_power_seed(seed + i)
    mu <- generate_power_signal(N_CAL, "updown", jump)
    y <- mu + rnorm(N_CAL)
    fit <- SVP(y, q, subtests = "both_lengthonly")
    truth <- normalise_boundaries(which(diff(mu) != 0), N_CAL)
    est <- normalise_boundaries(fit$changepoints, N_CAL)
    unname(paper_metrics(truth, est, TOLERANCE)["F1"])
  }, workers = workers)
  mean(unlist(scores))
}

cat("=== Bracket search ===\n")
q_lo <- 0.1; q_hi <- 3.0
reps <- 2000L; seed <- 740000L
fpr_lo <- empirical_fpr(q_lo, reps, seed)
fpr_hi <- empirical_fpr(q_hi, reps, seed)
cat(sprintf("q=%.3f -> FPR=%.4f | q=%.3f -> FPR=%.4f\n", q_lo, fpr_lo, q_hi, fpr_hi))
while (fpr_hi > ALPHA_TARGET && q_hi < 20) {
  q_hi <- q_hi * 1.5
  fpr_hi <- empirical_fpr(q_hi, reps, seed)
  cat(sprintf("widening: q=%.3f -> FPR=%.4f\n", q_hi, fpr_hi))
}
while (fpr_lo < ALPHA_TARGET && q_lo > 1e-6) {
  q_lo <- q_lo / 2
  fpr_lo <- empirical_fpr(q_lo, reps, seed)
  cat(sprintf("narrowing: q=%.3f -> FPR=%.4f\n", q_lo, fpr_lo))
}

cat("\n=== Bisecting q for a true 1% empirical FPR ===\n")
best <- list(q = q_hi, fpr = fpr_hi)
for (iter in 1:16) {
  q_mid <- (q_lo + q_hi) / 2
  fpr_mid <- empirical_fpr(q_mid, reps, seed)
  cat(sprintf("iter %2d: q=%.4f -> FPR=%.4f\n", iter, q_mid, fpr_mid))
  if (abs(fpr_mid - ALPHA_TARGET) < abs(best$fpr - ALPHA_TARGET)) best <- list(q = q_mid, fpr = fpr_mid)
  if (abs(fpr_mid - ALPHA_TARGET) < 0.0008) break
  if (fpr_mid > ALPHA_TARGET) q_lo <- q_mid else q_hi <- q_mid
}
q_star <- best$q
cat(sprintf("\nSelected q = %.4f (bisection FPR estimate = %.4f)\n", q_star, best$fpr))

cat("\n=== High-precision final validation (8000 reps) ===\n")
final_fpr <- empirical_fpr(q_star, 8000L, 960000L)
cat(sprintf("Final empirical FPR at n=1000 (8000 reps): %.4f\n", final_fpr))

cat("\n=== F1 at jump=0.6 on updown (2000 reps) ===\n")
final_f1 <- mean_updown_f1(q_star, jump = 0.6, reps = 2000L, seed = 1010000L)
cat(sprintf("F1@0.6 = %.4f\n", final_f1))

result <- list(q = q_star, final_fpr = final_fpr, final_f1_at_06 = final_f1)
saveRDS(result, file.path("simulations", "power_gaussian",
                          "lengthonly_calibration.rds"))
print(result)
