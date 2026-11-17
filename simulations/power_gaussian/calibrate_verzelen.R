## Exact multiscale critical values; gamma is a direct threshold.
source(file.path("simulations", "power_gaussian", "calibrate_modes.R"))

run_gaussian_verzelen_calibration <- function(
    n = 1000L, calibration_reps = 20000L, validation_reps = 10000L,
    calibration_seed = 810000L, validation_seed = 910000L,
    target = 0.01, workers = gaussian_calibration_workers(),
    root = GAUSSIAN_CALIBRATION_ROOT) {
  if (length(intersect(calibration_seed + seq_len(calibration_reps),
                       validation_seed + seq_len(validation_reps))))
    stop("Calibration and validation seed ranges must be disjoint")
  critical_values <- function(reps, seed) {
    values <- calibration_lapply(seq_len(reps), function(i) {
      y <- simulate_gaussian_null(n, seed + i)
      svpChange:::.verzelen_statistic_cpp(y)
    }, workers)
    if (!all(vapply(values, function(x)
        is.numeric(x) && length(x) == 1L && is.finite(x), logical(1))))
      stop("At least one critical-value computation failed")
    unlist(values, use.names = FALSE)
  }
  message("Computing ", calibration_reps, " Verzelen calibration maxima")
  calibration <- critical_values(calibration_reps, calibration_seed)
  center <- as.numeric(stats::quantile(calibration, 1 - target, type = 1))
  grid <- seq(floor((center - 0.05) * 1000),
              ceiling((center + 0.05) * 1000)) / 1000
  summarize <- function(values, gammas, stage, seed) {
    do.call(rbind, lapply(gammas, function(gamma) {
      failures <- sum(values > gamma)
      interval <- wilson_interval(failures, length(values))
      data.frame(gamma = gamma, false_detection = failures / length(values),
                 lower = interval[["lower"]], upper = interval[["upper"]],
                 reps = length(values), n = n, seed = seed, target_fpr = target,
                 stage = stage, stringsAsFactors = FALSE)
    }))
  }
  output <- summarize(calibration, grid, "calibration", calibration_seed)
  gap <- abs(output$false_detection - target)
  gamma <- max(output$gamma[gap == min(gap)])
  output$selected <- output$gamma == gamma
  message("Selected gamma = ", gamma, "; computing ", validation_reps,
          " independent validation maxima")
  validation <- critical_values(validation_reps, validation_seed)
  checked <- summarize(validation, gamma, "validation", validation_seed)
  checked$selected <- TRUE
  output <- rbind(output, checked)
  output$sigma2 <- 1
  output$formula <- "max_all_triples C^2/(4*sigma2)-2*log(N*(c-a)/((b-a)*(c-b)))"
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(calibration = calibration, validation = validation,
               gamma = gamma, n = n, calibration_seed = calibration_seed,
               validation_seed = validation_seed, target_fpr = target,
               formula = output$formula[1L]),
          file.path(root, "gaussian_verzelen_calibration_raw.rds"))
  utils::write.csv(output, file.path(root, "gaussian_verzelen_calibration.csv"),
                   row.names = FALSE)
  output[output$selected, ]
}

if (identical(tolower(Sys.getenv("SVP_RUN_VERZELEN_CALIBRATION")), "true"))
  print(run_gaussian_verzelen_calibration())
