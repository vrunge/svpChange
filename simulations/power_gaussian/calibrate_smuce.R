## Calibrate the SMUCE validity threshold to 1% null false detection.
library(svpChange)
library(Rcpp)

source(file.path("simulations", "power_gaussian", "calibrate_modes.R"))

GAUSSIAN_SMUCE_CALIBRATION_FILE <- file.path(
  "simulations", "power_gaussian", "gaussian_smuce_calibration.csv"
)
GAUSSIAN_SMUCE_CALIBRATION_REPS <- 20000L
GAUSSIAN_SMUCE_VALIDATION_REPS <- 10000L
GAUSSIAN_SMUCE_CALIBRATION_SEED <- 11100000L
GAUSSIAN_SMUCE_VALIDATION_SEED <- 12100000L

compile_smuce_critical_values <- function() {
  Rcpp::sourceCpp(file.path(
    "simulations", "power_gaussian", "smuce_critical_values.cpp"
  ))
}

run_smuce_critical_sample <- function(n, reps, seed, workers) {
  rows <- calibration_lapply(seq_len(reps), function(i) {
    y <- simulate_gaussian_null(n, seed + i)
    smuce_critical_value_cpp(y, sigma2 = 1)
  }, workers)
  values <- unlist(rows, use.names = FALSE)
  if (length(values) != reps || !all(is.finite(values)) || any(values < 0))
    stop("At least one SMUCE critical-value computation failed")
  values
}

run_gaussian_smuce_calibration <- function(
    n = GAUSSIAN_CALIBRATION_N,
    calibration_reps = GAUSSIAN_SMUCE_CALIBRATION_REPS,
    validation_reps = GAUSSIAN_SMUCE_VALIDATION_REPS,
    calibration_seed = GAUSSIAN_SMUCE_CALIBRATION_SEED,
    validation_seed = GAUSSIAN_SMUCE_VALIDATION_SEED,
    target = GAUSSIAN_FALSE_DETECTION_TARGET,
    workers = gaussian_calibration_workers(),
    root = dirname(GAUSSIAN_SMUCE_CALIBRATION_FILE)) {
  if (length(intersect(calibration_seed + seq_len(calibration_reps),
                       validation_seed + seq_len(validation_reps))))
    stop("Calibration and validation seed ranges must be disjoint")
  compile_smuce_critical_values()
  message("Computing ", calibration_reps, " SMUCE calibration thresholds")
  calibration <- run_smuce_critical_sample(
    n, calibration_reps, calibration_seed, workers
  )
  gamma <- as.numeric(stats::quantile(
    calibration, probs = 1 - target, type = 1, names = FALSE
  ))
  calibration_false <- sum(calibration > gamma)
  calibration_interval <- wilson_interval(calibration_false, calibration_reps)

  message("Selected gamma = ", format(gamma, digits = 8),
          "; computing ", validation_reps, " independent validation thresholds")
  validation <- run_smuce_critical_sample(
    n, validation_reps, validation_seed, workers
  )
  validation_false <- sum(validation > gamma)
  validation_interval <- wilson_interval(validation_false, validation_reps)
  output <- data.frame(
    gamma = gamma,
    target_fpr = target,
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
    n = n,
    sigma2 = 1,
    statistic = "minimum gamma making the full-series SMUCE constraint feasible",
    stringsAsFactors = FALSE
  )
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(output, GAUSSIAN_SMUCE_CALIBRATION_FILE, row.names = FALSE)
  saveRDS(list(calibration = calibration, validation = validation,
               result = output),
          sub("\\.csv$", ".rds", GAUSSIAN_SMUCE_CALIBRATION_FILE))
  output
}

verify_smuce_critical_values <- function(seed = 13100000L, reps = 100L) {
  compile_smuce_critical_values()
  for (i in seq_len(reps)) {
    y <- simulate_gaussian_null(80L, seed + i)
    critical <- smuce_critical_value_cpp(y, sigma2 = 1)
    above <- svp_smuce_cpp(y, q = critical + 1e-7, sigma2 = 1)
    if (length(above) != 1L)
      stop("SMUCE critical value failed its one-segment upper check")
    if (critical > 1e-6) {
      below <- svp_smuce_cpp(y, q = critical - 1e-7, sigma2 = 1)
      if (length(below) == 1L)
        stop("SMUCE critical value failed its one-segment lower check")
    }
  }
  invisible(TRUE)
}

if (identical(tolower(Sys.getenv("SVP_RUN_SMUCE_CALIBRATION")), "true")) {
  verify_smuce_critical_values()
  print(run_gaussian_smuce_calibration())
}
