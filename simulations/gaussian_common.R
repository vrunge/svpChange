## Shared names and calibrated constants for the Gaussian paper simulations.

GAUSSIAN_CALIBRATION_FILE <- file.path(
  "simulations", "power_gaussian", "gaussian_mode_calibration.csv"
)

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
