## Calibration for "both_empirical": a time-varying threshold gamma(t) for
## subtests="both", built directly from simulation rather than any
## parametric or asymptotic formula.
##
## gaussian_both_trajectory_cpp returns, for each prefix length b = 2..n,
## M_b = max_{0<=u<tau<b} Z(u,tau,b)^2: the RIGHT ENDPOINT IS FIXED AT b,
## only the start u and internal split tau vary. This is F_left(y_{0..b})
## in the F_right/F_left/F_both notation of the SVP chapter (Runge et al.,
## "Exact Multiple Change-Point Detection via Smallest Valid Partitioning"):
## max over every start, right endpoint fixed.
##
## The statistic that must actually be controlled is F_b = F_both(y_{0..b})
## = max over EVERY sub-segment of y_{0..b} (both endpoints free), because
## "both" pruning's discard-on-failure AND inclusive-left-pruning are only
## jointly exact when the validity test is F_both itself (gamma- AND
## gamma^{-1}-stable). Since F_both(y_{0..b}) = max_{b'<=b} F_left(y_{0..b'})
## = cummax(M)_b, F_b is obtained from the raw M trajectory by a running
## cummax over b -- exactly matching how gaussian_critical_values_cpp's own
## "$both" (the flat-threshold calibration) scans every start AND every
## intermediate endpoint, not just those ending at the final n.
##
## For b = 2..n, E_b/V_b are estimated by simulation on the cummax'd F_b,
## and G_b = (F_b - E_b)/sqrt(V_b) standardizes each. c is calibrated from
## the empirical 99th percentile of max_{b=2..n} G_b, using the JOINT
## (correlated) distribution of one full running trajectory per replicate --
## not a per-length or Bonferroni approximation. The test is then
## F_b > E_b + c*sqrt(V_b), applied via "both" pruning: every simultaneously
## alive candidate at time t=b compares its own max_cost() against the same
## shared threshold gamma(b), exactly like a flat "both" threshold but
## varying over time instead of being fixed once at t=n.
source(file.path("simulations", "power_gaussian", "calibrate_modes.R"))

compile_gaussian_both_trajectory <- function() {
  include <- normalizePath("src", mustWork = TRUE)
  old_flags <- Sys.getenv("PKG_CPPFLAGS", unset = NA_character_)
  on.exit({
    if (is.na(old_flags)) Sys.unsetenv("PKG_CPPFLAGS") else
      Sys.setenv(PKG_CPPFLAGS = old_flags)
  }, add = TRUE)
  Sys.setenv(PKG_CPPFLAGS = paste("-I", shQuote(include)))
  Rcpp::sourceCpp(file.path("simulations", "power_gaussian",
                            "gaussian_critical_values.cpp"))
}

N_CAL <- 1000L
REPS <- 20000L
SEED <- 5310000L
ALPHA <- 0.01

compile_gaussian_both_trajectory()
workers <- gaussian_calibration_workers()
if (workers < 2L) workers <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
cat("workers:", workers, "\n")

cat("Simulating", REPS, "running trajectories (n =", N_CAL, ") ...\n")
t0 <- Sys.time()
trajectories <- calibration_lapply(seq_len(REPS), function(i) {
  y <- simulate_gaussian_null(N_CAL, SEED + i)
  gaussian_both_trajectory_cpp(y) # length n-1, M_2..M_n = F_left(y_{0..b})
}, workers)
mat <- do.call(rbind, trajectories) # REPS x (n-1), raw F_left trajectories
cat("done in", as.numeric(Sys.time() - t0, units = "secs"), "s; dim:",
    paste(dim(mat), collapse = " x "), "\n")

saveRDS(mat, file.path("simulations", "power_gaussian",
                       "empirical_both_trajectories_Fleft_raw.rds"))

## F_both(y_{0..b}) = max_{b'<=b} F_left(y_{0..b'}): running cummax over b,
## row by row (each row is one replicate's trajectory).
mat <- t(apply(mat, 1, cummax))

saveRDS(mat, file.path("simulations", "power_gaussian",
                       "empirical_both_trajectories.rds"))

## Per-length moments (columns = b = 2..n)
E_b <- colMeans(mat)
V_b <- apply(mat, 2, var)
b_values <- seq(2, N_CAL)

cat("Standardizing and computing the joint running max ...\n")
G <- sweep(mat, 2, E_b, "-")
G <- sweep(G, 2, sqrt(V_b), "/")
running_max <- apply(G, 1, max)

c_value <- as.numeric(quantile(running_max, 1 - ALPHA, type = 1))
cat(sprintf("c (alpha=%.2f, %d reps) = %.6f\n", ALPHA, REPS, c_value))

calibration_fpr <- mean(running_max > c_value)
cat(sprintf("calibration-sample exceedance rate at c: %.4f (target %.2f)\n",
            calibration_fpr, ALPHA))

result <- list(b = b_values, E_b = E_b, V_b = V_b, c = c_value,
              n = N_CAL, reps = REPS, seed = SEED, alpha = ALPHA)
saveRDS(result, file.path("simulations", "power_gaussian",
                          "empirical_both_calibration.rds"))
utils::write.csv(
  data.frame(b = b_values, E_b = E_b, V_b = V_b),
  file.path("simulations", "power_gaussian", "empirical_both_calibration.csv"),
  row.names = FALSE
)
cat("Saved calibration table and c.\n")
