## Calibration for "both_correction": the multiscale "both" validity test of
## simulations/consistency/SVP_consistency_sharp.tex, Section 6
## ("Multiscale validity with both subsegment tests"). gamma (= q in the
## paper) is a direct threshold, not a log(n) multiplier. Calibration follows
## the paper's own recipe (Lemma "Simultaneous CUSUM control"): simulate the
## null maximum M_n = T_n(sigma*eps;0,n) and take its empirical 1-alpha
## quantile.
source(file.path("simulations", "power_gaussian", "calibrate_modes.R"))

GAUSSIAN_CONSISTENCY_SHARP_CALIBRATION_FILE <- file.path(
  "simulations", "power_gaussian", "gaussian_consistency_sharp_calibration.csv"
)

run_gaussian_consistency_sharp_calibration <- function(
    n = 1000L, calibration_reps = 20000L, validation_reps = 10000L,
    calibration_seed = 1310000L, validation_seed = 1410000L,
    target = 0.01, workers = gaussian_calibration_workers(),
    root = GAUSSIAN_CALIBRATION_ROOT) {
  if (length(intersect(calibration_seed + seq_len(calibration_reps),
                       validation_seed + seq_len(validation_reps))))
    stop("Calibration and validation seed ranges must be disjoint")
  critical_values <- function(reps, seed) {
    values <- calibration_lapply(seq_len(reps), function(i) {
      y <- simulate_gaussian_null(n, seed + i)
      svpChange:::.consistency_sharp_statistic_cpp(y)
    }, workers)
    if (!all(vapply(values, function(x)
        is.numeric(x) && length(x) == 1L && is.finite(x), logical(1))))
      stop("At least one critical-value computation failed")
    unlist(values, use.names = FALSE)
  }
  message("Computing ", calibration_reps, " both_correction calibration maxima")
  calibration <- critical_values(calibration_reps, calibration_seed)
  q <- as.numeric(stats::quantile(calibration, 1 - target, type = 1))
  calibration_false <- sum(calibration > q)
  calibration_interval <- wilson_interval(calibration_false, calibration_reps)

  message("Selected q = ", format(q, digits = 8), "; computing ",
          validation_reps, " independent validation maxima")
  validation <- critical_values(validation_reps, validation_seed)
  validation_false <- sum(validation > q)
  validation_interval <- wilson_interval(validation_false, validation_reps)

  output <- data.frame(
    q = q, target_fpr = target,
    calibration_false_detections = calibration_false,
    calibration_fpr = calibration_false / calibration_reps,
    calibration_lower = calibration_interval[["lower"]],
    calibration_upper = calibration_interval[["upper"]],
    calibration_reps = calibration_reps, calibration_seed = calibration_seed,
    validation_false_detections = validation_false,
    validation_fpr = validation_false / validation_reps,
    validation_lower = validation_interval[["lower"]],
    validation_upper = validation_interval[["upper"]],
    validation_reps = validation_reps, validation_seed = validation_seed,
    n = n, sigma2 = 1,
    formula = "max_{a<u<b} |Z_u(a,b)| - sqrt(2*log(n*(b-a)/((u-a)*(b-u))))",
    stringsAsFactors = FALSE
  )
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(calibration = calibration, validation = validation, result = output),
          sub("\\.csv$", ".rds", GAUSSIAN_CONSISTENCY_SHARP_CALIBRATION_FILE))
  utils::write.csv(output, GAUSSIAN_CONSISTENCY_SHARP_CALIBRATION_FILE, row.names = FALSE)
  output
}

read_gaussian_consistency_sharp_calibration <- function(
    path = GAUSSIAN_CONSISTENCY_SHARP_CALIBRATION_FILE, n = 1000L) {
  values <- utils::read.csv(path, stringsAsFactors = FALSE)
  if (nrow(values) != 1L || !is.finite(values$q) || values$n != n)
    stop("both_correction calibration must contain one threshold for n = ", n)
  values$q
}

if (identical(tolower(Sys.getenv("SVP_RUN_CONSISTENCY_SHARP_CALIBRATION")), "true"))
  print(run_gaussian_consistency_sharp_calibration())
