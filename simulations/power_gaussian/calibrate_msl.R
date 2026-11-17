## Calibrate exact SVP both MS-l validity to 1% no-change detection.
library(svpChange)
library(Rcpp)

source(file.path("simulations", "power_common.R"))
source(file.path("simulations", "gaussian_common.R"))

GAUSSIAN_MSL_CALIBRATION_FILE <- file.path(
  "simulations", "power_gaussian", "gaussian_msl_calibration.csv"
)

compile_msl_critical_values <- function() {
  include <- normalizePath("src", mustWork = TRUE)
  old_flags <- Sys.getenv("PKG_CPPFLAGS", unset = NA_character_)
  on.exit({
    if (is.na(old_flags)) Sys.unsetenv("PKG_CPPFLAGS") else
      Sys.setenv(PKG_CPPFLAGS = old_flags)
  }, add = TRUE)
  Sys.setenv(PKG_CPPFLAGS = paste("-I", shQuote(include)))
  Rcpp::sourceCpp(file.path(
    "simulations", "power_gaussian", "msl_critical_values.cpp"
  ))
}

msl_wilson_interval <- function(successes, trials, level = 0.95) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  p <- successes / trials
  denominator <- 1 + z^2 / trials
  center <- (p + z^2 / (2 * trials)) / denominator
  half_width <- z * sqrt(p * (1 - p) / trials +
                           z^2 / (4 * trials^2)) / denominator
  c(lower = max(0, center - half_width),
    upper = min(1, center + half_width))
}

run_msl_critical_sample <- function(n, reps, seed, sigma2, workers) {
  values <- power_lapply(seq_len(reps), function(i) {
    set_power_seed(seed + i)
    gaussian_msl_critical_value_cpp(stats::rnorm(n, sd = sqrt(sigma2)),
                                    sigma2)
  }, workers = workers)
  values <- unlist(values, use.names = FALSE)
  if (length(values) != reps || any(!is.finite(values)))
    stop("At least one MS-l critical-value computation failed")
  values
}

verify_msl_critical_values <- function(seed = 61000000L, reps = 30L) {
  compile_msl_critical_values()
  for (i in seq_len(reps)) {
    set_power_seed(seed + i)
    sigma2 <- if (i %% 2L) 1 else 2.25
    y <- stats::rnorm(40L, sd = sqrt(sigma2))
    critical <- gaussian_msl_critical_value_cpp(y, sigma2)
    above <- SVP(y, critical + 1e-7, subtests = "both_msl",
                 sigma2 = sigma2)$changepoints
    below <- SVP(y, critical - 1e-7, subtests = "both_msl",
                 sigma2 = sigma2)$changepoints
    if (length(above) != 1L || length(below) == 1L)
      stop("MS-l critical value disagrees with the SVP one-segment event")
  }
  invisible(TRUE)
}

run_gaussian_msl_calibration <- function(
    n = 1000L, calibration_reps = 20000L, validation_reps = 10000L,
    calibration_seed = 41000000L, validation_seed = 51000000L,
    sigma2 = 1, target = 0.01, workers = power_default_workers()) {
  if (length(n) != 1L || !is.finite(n) || n < 2L || n != as.integer(n))
    stop("n must be an integer of at least two")
  if (length(calibration_reps) != 1L || calibration_reps < 1L ||
      length(validation_reps) != 1L || validation_reps < 1L)
    stop("replication counts must be positive")
  if (length(sigma2) != 1L || !is.finite(sigma2) || sigma2 <= 0)
    stop("sigma2 must be finite and positive")
  if (length(target) != 1L || !is.finite(target) ||
      target <= 0 || target >= 1)
    stop("target must lie in (0, 1)")
  if (calibration_seed + calibration_reps >= validation_seed &&
      validation_seed + validation_reps >= calibration_seed)
    stop("calibration and validation seed ranges must be disjoint")

  compile_msl_critical_values()
  message("Computing ", calibration_reps, " MS-l calibration maxima")
  calibration <- run_msl_critical_sample(
    n, calibration_reps, calibration_seed, sigma2, workers
  )
  q <- as.numeric(stats::quantile(
    calibration, probs = 1 - target, type = 1, names = FALSE
  ))
  calibration_false <- sum(calibration > q)
  calibration_interval <- msl_wilson_interval(calibration_false,
                                               calibration_reps)

  message("Selected q = ", format(q, digits = 8),
          "; computing ", validation_reps, " independent maxima")
  validation <- run_msl_critical_sample(
    n, validation_reps, validation_seed, sigma2, workers
  )
  validation_false <- sum(validation > q)
  validation_interval <- msl_wilson_interval(validation_false,
                                              validation_reps)
  output <- data.frame(
    q = q, target_fpr = target,
    calibration_false_detections = calibration_false,
    calibration_fpr = calibration_false / calibration_reps,
    calibration_lower = calibration_interval[["lower"]],
    calibration_upper = calibration_interval[["upper"]],
    calibration_reps = calibration_reps,
    calibration_seed = calibration_seed,
    validation_false_detections = validation_false,
    validation_fpr = validation_false / validation_reps,
    validation_lower = validation_interval[["lower"]],
    validation_upper = validation_interval[["upper"]],
    validation_reps = validation_reps,
    validation_seed = validation_seed,
    n = n, sigma2 = sigma2,
    statistic = "max over intervals and splits of abs(Z)/sigma - sqrt(2 log(e length / shorter_side))",
    stringsAsFactors = FALSE
  )
  utils::write.csv(output, GAUSSIAN_MSL_CALIBRATION_FILE, row.names = FALSE)
  saveRDS(list(calibration = calibration, validation = validation,
               result = output),
          sub("\\.csv$", ".rds", GAUSSIAN_MSL_CALIBRATION_FILE))
  output
}

if (identical(tolower(Sys.getenv("SVP_CALIBRATE_GAUSSIAN_MSL")), "true")) {
  verify_msl_critical_values()
  print(run_gaussian_msl_calibration())
}
