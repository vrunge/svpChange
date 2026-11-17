## Shared names and calibrated constants for the Gaussian paper simulations.

GAUSSIAN_CALIBRATION_FILE <- file.path(
  "simulations", "power_gaussian", "gaussian_mode_calibration.csv"
)
GAUSSIAN_PELT_CALIBRATION_FILE <- file.path(
  "simulations", "power_gaussian", "gaussian_pelt_calibration.csv"
)
GAUSSIAN_VERZELEN_CALIBRATION_FILE <- file.path(
  "simulations", "power_gaussian", "gaussian_verzelen_calibration.csv"
)
GAUSSIAN_SMUCE_CALIBRATION_FILE <- file.path(
  "simulations", "power_gaussian", "gaussian_smuce_calibration.csv"
)
GAUSSIAN_SMUCE_UP_TUNING_FILE <- file.path(
  "simulations", "power_gaussian", "gaussian_smuce_up_tuning.csv"
)
GAUSSIAN_DYADIC_CALIBRATION_FILE <- file.path(
  "simulations", "power_gaussian", "gaussian_dyadic_calibration.csv"
)
GAUSSIAN_MSL_CALIBRATION_FILE <- file.path(
  "simulations", "power_gaussian", "gaussian_msl_calibration.csv"
)
GAUSSIAN_CONSISTENCY_SHARP_CALIBRATION_FILE <- file.path(
  "simulations", "power_gaussian", "gaussian_consistency_sharp_calibration.csv"
)

read_gaussian_consistency_sharp_calibration <- function(
    path = GAUSSIAN_CONSISTENCY_SHARP_CALIBRATION_FILE, n = 1000L) {
  values <- utils::read.csv(path, stringsAsFactors = FALSE)
  if (nrow(values) != 1L || !is.finite(values$q) || values$n != n)
    stop("both_correction calibration must contain one threshold for n = ", n)
  values$q
}

## Darling-Erdos-type sharp null threshold for subtests = "none"
## (SVP_consistency_sharp.tex, Theorem 1 and the "Fixed error probability"
## remark): d_n = log(log(n)) + log(log(log(n)))/2 - log(pi)/2, and the
## asymptotic alpha-level threshold is d_n + log(2) - log(-log(1 - alpha)) on
## the paper's cost scale f = (1/2) max Z^2. The package's raw statistic
## (max_cost(), no 1/2 factor) is compared against twice that value.
gaussian_sharp_none_threshold <- function(n, alpha = 0.01) {
  d_n <- log(log(n)) + 0.5 * log(log(log(n))) - 0.5 * log(pi)
  gamma_paper <- d_n + log(2) - log(-log(1 - alpha))
  2 * gamma_paper
}

read_gaussian_msl_calibration <- function(
    path = GAUSSIAN_MSL_CALIBRATION_FILE, n = 1000L) {
  values <- utils::read.csv(path, stringsAsFactors = FALSE)
  if (nrow(values) != 1L || !is.finite(values$q) ||
      values$n != n || values$sigma2 != 1)
    stop("MS-l calibration must contain one threshold for n = ", n,
         " and sigma2 = 1")
  values$q
}

read_gaussian_verzelen_calibration <- function(
    path = GAUSSIAN_VERZELEN_CALIBRATION_FILE, n = 1000L) {
  values <- utils::read.csv(path, stringsAsFactors = FALSE)
  selected <- values[values$stage == "calibration" & values$selected, , drop = FALSE]
  if (nrow(selected) != 1L || !is.finite(selected$gamma) || selected$n != n)
    stop("Verzelen calibration must contain one selected threshold for n = ", n)
  selected$gamma
}

read_gaussian_smuce_calibration <- function(
    path = GAUSSIAN_SMUCE_CALIBRATION_FILE, n = 1000L) {
  values <- utils::read.csv(path, stringsAsFactors = FALSE)
  if (nrow(values) != 1L || !is.finite(values$gamma) || values$n != n)
    stop("SMUCE calibration must contain one threshold for n = ", n)
  values$gamma
}

read_gaussian_smuce_up_tuning <- function(
    path = GAUSSIAN_SMUCE_UP_TUNING_FILE, n = 1000L) {
  values <- utils::read.csv(path, stringsAsFactors = FALSE)
  selected <- values[values$selected, , drop = FALSE]
  if (nrow(selected) != 1L || !is.finite(selected$gamma) || selected$n != n)
    stop("SMUCE up-tuned gamma must contain one selected value for n = ", n)
  selected$gamma
}

read_gaussian_dyadic_calibration <- function(
    path = GAUSSIAN_DYADIC_CALIBRATION_FILE, n = 1000L) {
  values <- utils::read.csv(path, stringsAsFactors = FALSE)
  selected <- values[values$selected, , drop = FALSE]
  selected <- selected[order(selected$band), , drop = FALSE]
  expected <- max(1L, floor(log2(n / 2)) + 1L)
  if (nrow(selected) != expected || any(!is.finite(selected$gamma)) ||
      any(selected$gamma <= 0) || any(selected$n != n) ||
      !identical(as.integer(selected$band), seq_len(expected) - 1L)) {
    stop("Dyadic calibration must contain one gamma per band for n = ", n)
  }
  selected$gamma
}

GAUSSIAN_UNCALIBRATED_RIGHT_CONSTANT <- 2

gaussian_power_method_labels <- function(constants, pelt_constant) {
  c(paste0("PELT / c = ", formatC(pelt_constant, digits = 3L, format = "f")),
    paste0("SVP none / c = ", formatC(constants[["none"]], digits = 3L, format = "f")),
    paste0("SVP right uncalibrated / c = ",
           formatC(GAUSSIAN_UNCALIBRATED_RIGHT_CONSTANT, digits = 3L, format = "f")),
    paste0("SVP right / c = ", formatC(constants[["right"]], digits = 3L, format = "f")),
    paste0("SVP both / c = ", formatC(constants[["both"]], digits = 3L, format = "f")))
}

validate_gaussian_power_results <- function(results) {
  expected <- gaussian_power_method_labels(
    read_gaussian_calibration(), read_gaussian_pelt_calibration())
  if (!setequal(unique(as.character(results$algorithm)), expected))
    stop("Stale Gaussian results: rerun the five-method power study ",
         "(PELT, SVP none/right-uncalibrated/right/both).")
  expected
}

gaussian_method_label <- function(mode, constant) {
  paste0(
    "SVP Gaussian / ", mode, " / c = ",
    formatC(constant, digits = 3L, format = "f")
  )
}

read_gaussian_calibration <- function(
    path = GAUSSIAN_CALIBRATION_FILE,
    required_modes = c("none", "right", "both")) {
  if (!file.exists(path)) {
    stop(
      "Gaussian calibration file not found: ", path,
      ". Run simulations/power_gaussian/calibrate_modes.R first."
    )
  }
  calibration <- utils::read.csv(path, stringsAsFactors = FALSE)
  required_columns <- c("mode", "constant")
  if (!all(required_columns %in% names(calibration))) {
    stop("Gaussian calibration must contain mode and constant columns")
  }
  missing_modes <- setdiff(required_modes, calibration$mode)
  if (length(missing_modes)) {
    stop("Missing Gaussian calibration for: ", paste(missing_modes, collapse = ", "))
  }
  calibration <- calibration[match(required_modes, calibration$mode), , drop = FALSE]
  stats::setNames(calibration$constant, calibration$mode)
}

read_gaussian_pelt_calibration <- function(
    path = GAUSSIAN_PELT_CALIBRATION_FILE) {
  if (!file.exists(path)) {
    stop("Gaussian PELT calibration file not found: ", path)
  }
  values <- utils::read.csv(path, stringsAsFactors = FALSE)
  selected <- values[
    values$stage == "calibration" & values$selected,
    "constant", drop = TRUE
  ]
  if (length(selected) != 1L) {
    stop("Gaussian PELT calibration must contain one selected constant")
  }
  as.numeric(selected)
}
