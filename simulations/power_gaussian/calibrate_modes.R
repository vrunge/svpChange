## Calibrate Gaussian SVP constants to the PELT null false-positive rate.

library(svpChange)
library(changepoint)
library(Rcpp)

source(file.path("simulations", "gaussian_common.R"))

GAUSSIAN_CALIBRATION_ROOT <- file.path("simulations", "power_gaussian")
GAUSSIAN_CALIBRATION_N <- 1000L
GAUSSIAN_CALIBRATION_REPS <- 20000L
GAUSSIAN_VALIDATION_REPS <- 10000L
GAUSSIAN_CALIBRATION_SEED <- 810000L
GAUSSIAN_VALIDATION_SEED <- 820000L

compile_gaussian_critical_values <- function() {
  include <- normalizePath("src", mustWork = TRUE)
  old_flags <- Sys.getenv("PKG_CPPFLAGS", unset = NA_character_)
  on.exit({
    if (is.na(old_flags)) Sys.unsetenv("PKG_CPPFLAGS") else
      Sys.setenv(PKG_CPPFLAGS = old_flags)
  }, add = TRUE)
  Sys.setenv(PKG_CPPFLAGS = paste("-I", shQuote(include)))
  Rcpp::sourceCpp(file.path(
    GAUSSIAN_CALIBRATION_ROOT, "gaussian_critical_values.cpp"
  ))
}

gaussian_calibration_workers <- function() {
  requested <- suppressWarnings(as.integer(Sys.getenv(
    "SVP_GAUSSIAN_CALIBRATION_WORKERS", "1"
  )))
  if (length(requested) != 1L || is.na(requested) || requested < 1L) 1L else requested
}

calibration_lapply <- function(X, FUN, workers) {
  if (workers == 1L || .Platform$OS.type == "windows") return(lapply(X, FUN))
  parallel::mclapply(X, FUN, mc.cores = workers, mc.preschedule = TRUE,
                     mc.set.seed = FALSE)
}

simulate_gaussian_null <- function(n, seed) {
  set.seed(seed, kind = "Mersenne-Twister",
           normal.kind = "Inversion", sample.kind = "Rejection")
  ts_generator(
    chpts = n, parameters = 0, sd_noise = 1, type = "gauss"
  )
}

pelt_has_false_positive <- function(y) {
  fit <- changepoint::cpt.mean(
    y, method = "PELT", penalty = "Manual",
    pen.value = 2 * log(length(y))
  )
  length(changepoint::cpts(fit)) > 0L
}

run_critical_sample <- function(n, reps, seed, workers) {
  rows <- calibration_lapply(seq_len(reps), function(i) {
    y <- simulate_gaussian_null(n, seed + i)
    critical <- gaussian_critical_values_cpp(y)
    c(
      rep = i,
      pelt_false_positive = as.numeric(pelt_has_false_positive(y)),
      critical
    )
  }, workers)
  values <- as.data.frame(do.call(rbind, rows))
  values$rep <- as.integer(values$rep)
  values
}

wilson_interval <- function(successes, trials, level = 0.95) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  p <- successes / trials
  denominator <- 1 + z^2 / trials
  center <- (p + z^2 / (2 * trials)) / denominator
  half_width <- z * sqrt(p * (1 - p) / trials + z^2 / (4 * trials^2)) /
    denominator
  c(lower = max(0, center - half_width), upper = min(1, center + half_width))
}

select_mode_constants <- function(values, n, step = 0.001) {
  target <- mean(values$pelt_false_positive)
  modes <- c("none", "right", "both")
  do.call(rbind, lapply(modes, function(mode) {
    ratios <- values[[mode]] / log(n)
    centre <- stats::quantile(ratios, probs = 1 - target, type = 1,
                              names = FALSE)
    grid <- seq(
      max(0.001, floor((centre - 0.05) / step) * step),
      ceiling((centre + 0.05) / step) * step,
      by = step
    )
    fpr <- vapply(grid, function(constant) {
      mean(values[[mode]] >= constant * log(n))
    }, numeric(1))
    gap <- abs(fpr - target)
    chosen <- max(which(gap == min(gap)))
    data.frame(
      mode = mode,
      constant = grid[chosen],
      pelt_fpr_calibration = target,
      svp_fpr_calibration = fpr[chosen],
      absolute_gap_calibration = gap[chosen],
      stringsAsFactors = FALSE
    )
  }))
}

validate_mode_constants <- function(constants, values, n) {
  target <- mean(values$pelt_false_positive)
  target_ci <- wilson_interval(sum(values$pelt_false_positive), nrow(values))
  rows <- lapply(seq_len(nrow(constants)), function(i) {
    mode <- constants$mode[i]
    false_positive <- values[[mode]] >= constants$constant[i] * log(n)
    interval <- wilson_interval(sum(false_positive), length(false_positive))
    data.frame(
      mode = mode,
      constant = constants$constant[i],
      pelt_fpr_validation = target,
      pelt_fpr_lower = target_ci[["lower"]],
      pelt_fpr_upper = target_ci[["upper"]],
      svp_fpr_validation = mean(false_positive),
      svp_fpr_lower = interval[["lower"]],
      svp_fpr_upper = interval[["upper"]],
      absolute_gap_validation = abs(mean(false_positive) - target),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

verify_critical_values <- function(constants, seed = 830000L, reps = 100L) {
  for (i in seq_len(reps)) {
    y <- simulate_gaussian_null(80L, seed + i)
    critical <- gaussian_critical_values_cpp(y)
    for (mode in constants$mode) {
      for (gamma in c(0.95, 1.05) * critical[[mode]] + 1e-8) {
        fit <- SVP(y, gamma, "gaussian_mean", subtests = mode)
        expected_one_segment <- critical[[mode]] < gamma
        if (identical(length(fit$changepoints) == 1L, expected_one_segment)) next
        stop("Critical-value verification failed for mode ", mode)
      }
    }
  }
  invisible(TRUE)
}

run_gaussian_mode_calibration <- function(
    n = GAUSSIAN_CALIBRATION_N,
    calibration_reps = GAUSSIAN_CALIBRATION_REPS,
    validation_reps = GAUSSIAN_VALIDATION_REPS,
    workers = gaussian_calibration_workers()) {
  compile_gaussian_critical_values()
  calibration_values <- run_critical_sample(
    n, calibration_reps, GAUSSIAN_CALIBRATION_SEED, workers
  )
  constants <- select_mode_constants(calibration_values, n)
  validation_values <- run_critical_sample(
    n, validation_reps, GAUSSIAN_VALIDATION_SEED, workers
  )
  validation <- validate_mode_constants(constants, validation_values, n)
  output <- merge(constants, validation, by = c("mode", "constant"), sort = FALSE)
  output <- output[match(c("none", "right", "both"), output$mode), ]
  output$n <- n
  output$calibration_reps <- calibration_reps
  output$validation_reps <- validation_reps
  output$calibration_seed <- GAUSSIAN_CALIBRATION_SEED
  output$validation_seed <- GAUSSIAN_VALIDATION_SEED

  verify_critical_values(output)
  saveRDS(calibration_values, file.path(
    GAUSSIAN_CALIBRATION_ROOT, "gaussian_mode_calibration_raw.rds"
  ))
  saveRDS(validation_values, file.path(
    GAUSSIAN_CALIBRATION_ROOT, "gaussian_mode_validation_raw.rds"
  ))
  utils::write.csv(output, GAUSSIAN_CALIBRATION_FILE, row.names = FALSE)
  output
}

if (identical(tolower(Sys.getenv("SVP_RUN_GAUSSIAN_CALIBRATION")), "true")) {
  print(run_gaussian_mode_calibration())
}
