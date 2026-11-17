## Head-to-head F1 comparison: "all" interval system SMUCE (svp_smuce_cpp,
## every sub-window) versus "dyaLen" interval system SMUCE
## (svp_smuce_dyadic_cpp, every starting position but only dyadic lengths),
## both calibrated to a 1% null false-detection rate at n = 1000. Same
## scenarios, jump sizes, reps and tolerance as
## run_and_save_gaussian_four_methods(), so the F1 curves are directly
## comparable to the main four-methods figure.
library(svpChange)
library(dplyr)

source(file.path("simulations", "power_common.R"))
source(file.path("simulations", "gaussian_common.R"))
GAUSSIAN_ROOT <- file.path("simulations", "power_gaussian")

GAUSSIAN_SMUCE_DYADIC_CALIBRATION_FILE <- file.path(
  "simulations", "power_gaussian", "gaussian_smuce_dyadic_calibration.csv"
)

read_gaussian_smuce_dyadic_calibration <- function(
    path = GAUSSIAN_SMUCE_DYADIC_CALIBRATION_FILE, n = 1000L) {
  values <- utils::read.csv(path, stringsAsFactors = FALSE)
  if (nrow(values) != 1L || !is.finite(values$gamma) || values$n != n)
    stop("Dyadic SMUCE calibration must contain one threshold for n = ", n)
  values$gamma
}

smuce_dyadic_method_labels <- function(
    smuce_gamma = read_gaussian_smuce_calibration(),
    smuce_dyadic_gamma = read_gaussian_smuce_dyadic_calibration()) {
  c(paste0("SVP SMUCE all / gamma = ", formatC(smuce_gamma, digits = 3L, format = "f")),
    paste0("SVP SMUCE dyadic / gamma = ", formatC(smuce_dyadic_gamma, digits = 3L, format = "f")))
}

fit_smuce_dyadic_comparison <- function(
    y, smuce_gamma = read_gaussian_smuce_calibration(n = length(y)),
    smuce_dyadic_gamma = read_gaussian_smuce_dyadic_calibration(n = length(y))) {
  n <- length(y)
  methods <- list(
    normalise_boundaries(svp_smuce_cpp(y, q = smuce_gamma, sigma2 = 1), n),
    normalise_boundaries(svp_smuce_dyadic_cpp(y, q = smuce_dyadic_gamma, sigma2 = 1), n))
  names(methods) <- smuce_dyadic_method_labels(smuce_gamma, smuce_dyadic_gamma)
  methods
}

run_and_save_smuce_dyadic_comparison <- function(
    n = 1000L, jump_sizes = seq(0.1, 2, 0.1), reps = 100L,
    workers = power_default_workers(), seed = 7793001L, root = GAUSSIAN_ROOT) {
  labels <- smuce_dyadic_method_labels()
  results <- run_power_grid(
    n, jump_sizes, reps, simulate_noise = stats::rnorm,
    fit_methods = function(y, true_segments) fit_smuce_dyadic_comparison(y),
    tolerance = round(n * 0.0025), workers = workers, seed = seed)
  results <- set_algorithm_order(results, labels)
  saveRDS(results, file.path(root, "smuce_dyadic_comparison_results.rds"))
  utils::write.csv(dplyr::select(results, -changepoints),
                   file.path(root, "smuce_dyadic_comparison_results.csv"),
                   row.names = FALSE)
  invisible(results)
}

if (identical(tolower(Sys.getenv("SVP_RUN_SMUCE_DYADIC_COMPARISON")), "true"))
  run_and_save_smuce_dyadic_comparison()
