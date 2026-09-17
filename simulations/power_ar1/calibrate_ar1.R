## Calibration for the corrected Section 5.3 AR(1) power study.

library(svpChange)
library(changepoint)

source(file.path("simulations", "power_common.R"))

AR1_ROOT <- file.path("simulations", "power_ar1")
AR1_N <- 600L
AR1_RHO <- 0.8
AR1_TOLERANCE <- 5L
AR1_CALIBRATION_REPS <- 1000L
AR1_VALIDATION_REPS <- 1000L
AR1_CALIBRATION_SEED <- 930000L
AR1_VALIDATION_SEED <- 940000L
AR1_NULL_TARGET <- 0.99

ar1_innovation_variance <- function(rho = AR1_RHO, marginal_sd = 1) {
  marginal_sd^2 * (1 - rho^2)
}

simulate_ar1_noise <- function(n, rho = AR1_RHO, marginal_sd = 1) {
  ts_generator(
    chpts = n,
    parameters = 0,
    sd_noise = marginal_sd * sqrt(1 - rho^2),
    rho = rho,
    type = "gaussAR1"
  )
}

simulate_ar1_null <- function(n, rho, seed) {
  set_power_seed(seed)
  simulate_ar1_noise(n, rho)
}

ar1_null_stream <- function(n, rho, reps, seed) {
  lapply(seq_len(reps), function(i) simulate_ar1_null(n, rho, seed + i))
}

# DeCAFS 3.3.6 fixed the weight given to the first observation in its dynamic
# programme: it used precision `1 / (sdNu^2 (1 - phi^2))` for `mu_1` instead of
# the stationary AR(1) precision `(1 - phi^2) / sdNu^2` that the model implies.
# At rho = 0.8 that over-weights y[1] by a factor 7.7, which buys an almost free
# spurious change in the first few observations and inflates the penalty needed
# to reach the null target from about 2.75*log(n) to 7.75*log(n). The study is
# not comparable against earlier builds, so the version is checked up front.
DECAFS_MIN_VERSION <- "3.3.6"

decafs_version <- function() {
  as.character(utils::packageVersion("DeCAFS"))
}

require_decafs <- function() {
  if (!requireNamespace("DeCAFS", quietly = TRUE)) {
    stop("the optional package 'DeCAFS' is required for the AR(1) study")
  }
  if (utils::packageVersion("DeCAFS") < DECAFS_MIN_VERSION) {
    stop("DeCAFS ", DECAFS_MIN_VERSION, " or later is required; found ",
         utils::packageVersion("DeCAFS"))
  }
  invisible(TRUE)
}

wilson_interval <- function(successes, trials, level = 0.95) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  p <- successes / trials
  denominator <- 1 + z^2 / trials
  centre <- (p + z^2 / (2 * trials)) / denominator
  half_width <- z * sqrt(p * (1 - p) / trials + z^2 / (4 * trials^2)) /
    denominator
  c(lower = max(0, centre - half_width),
    upper = min(1, centre + half_width))
}

pelt_ar1_approximate_boundaries <- function(
    y, rho = AR1_RHO, penalty = 2.2 * log(length(y))) {
  n <- length(y)
  transformed <- y[-1L] - rho * y[-n]
  fit <- svpChange::PELT(transformed, penalty = penalty)
  internal <- head(fit$changepoints, -1L)
  normalise_boundaries(internal + 1L, n)
}

pelt_inflated_boundaries <- function(
    y, rho = AR1_RHO, penalty = 3 * log(length(y)) *
      (1 + rho) / (1 - rho)) {
  n <- length(y)
  fit <- changepoint::cpt.mean(
    y, method = "PELT", penalty = "Manual", pen.value = penalty
  )
  normalise_boundaries(changepoint::cpts(fit), n)
}

decafs_ar1_boundaries <- function(
    y, rho = AR1_RHO, penalty = 1 * log(length(y))) {
  n <- length(y)
  fit <- DeCAFS::DeCAFS(
    y,
    beta = penalty,
    modelParam = list(
      sdEta = 0,
      sdNu = sqrt(ar1_innovation_variance(rho)),
      phi = rho
    )
  )
  normalise_boundaries(fit$changepoints, n)
}

# SVP ranks the partitions that attain the smallest number of valid segments by
# their within-segment cost, so `cost` moves the boundaries without changing how
# many there are. "ar1" scores each segment by the sum of squared deviations of
# its interior innovations z[u] = y[u] - rho*y[u-1], skipping the innovation
# that straddles the boundary because its mean carries the jump. "gaussian" is
# the squared-error cost on the raw observations, which is what this study used
# before and which is misspecified under serial dependence.
svp_ar1focus_boundaries <- function(
    y, rho = AR1_RHO, constant = 1, cost = "ar1") {
  n <- length(y)
  normalise_boundaries(
    SVP(
      y, constant * log(n), "AR1Focus", subtests = "right",
      rho = rho, sigma2 = ar1_innovation_variance(rho), cost = cost
    )$changepoints,
    n
  )
}

ar1_fixed_pelt_boundaries <- function(y, method, rho = AR1_RHO) {
  switch(
    method,
    "PELT AR1 approximate" = pelt_ar1_approximate_boundaries(y, rho),
    "PELT inflated" = pelt_inflated_boundaries(y, rho),
    stop("Unknown fixed AR(1) PELT method: ", method)
  )
}

ar1_null_f1 <- function(boundaries, n, tolerance = AR1_TOLERANCE) {
  paper_metrics(c(n), boundaries, tolerance, terminal = n)[["F1"]]
}

fit_ar1_calibration_candidate <- function(
    y, method, candidate, n, rho) {
  switch(
    method,
    "DeCAFS AR1" = decafs_ar1_boundaries(
      y, rho, penalty = candidate * log(n)
    ),
    "PELT AR1 approximate" = pelt_ar1_approximate_boundaries(
      y, rho, penalty = candidate * log(n)
    ),
    "SVP AR1Focus" = svp_ar1focus_boundaries(
      y, rho, constant = candidate
    ),
    stop("Unknown AR(1) calibration method: ", method)
  )
}

evaluate_ar1_candidates <- function(
    candidates, method, null_series, n, rho, tolerance = AR1_TOLERANCE,
    workers = 1L) {
  score_one <- function(candidate) {
    scores <- vapply(null_series, function(y) {
      boundaries <- fit_ar1_calibration_candidate(
        y, method, candidate, n, rho
      )
      ar1_null_f1(boundaries, n, tolerance)
    }, numeric(1))
    interval <- wilson_interval(sum(scores == 1), length(scores))
    data.frame(
      method = method,
      candidate = candidate,
      null_f1 = mean(scores),
      null_f1_lower = interval[["lower"]],
      null_f1_upper = interval[["upper"]],
      recovered_one_segment = sum(scores == 1),
      mean_segments = mean(vapply(
        null_series,
        function(y) length(fit_ar1_calibration_candidate(
          y, method, candidate, n, rho
        )),
        integer(1)
      )),
      reps = length(null_series),
      n = n,
      rho = rho,
      tolerance = tolerance,
      decafs_version = decafs_version(),
      stringsAsFactors = FALSE
    )
  }
  rows <- power_lapply(candidates, score_one, workers = workers)
  dplyr::bind_rows(rows)
}

evaluate_fixed_ar1_method <- function(
    method, null_series, n, rho, seed, tolerance = AR1_TOLERANCE) {
  scores <- vapply(null_series, function(y) {
    ar1_null_f1(ar1_fixed_pelt_boundaries(y, method, rho), n, tolerance)
  }, numeric(1))
  interval <- wilson_interval(sum(scores == 1), length(scores))
  data.frame(
    method = method,
    candidate = if (method == "PELT AR1 approximate") 2.2 else 3,
    null_f1 = mean(scores),
    null_f1_lower = interval[["lower"]],
    null_f1_upper = interval[["upper"]],
    recovered_one_segment = sum(scores == 1),
    mean_segments = mean(vapply(
      null_series, function(y) length(ar1_fixed_pelt_boundaries(y, method, rho)),
      integer(1)
    )),
    reps = length(null_series),
    n = n,
    rho = rho,
    tolerance = tolerance,
    decafs_version = decafs_version(),
    seed = seed,
    stage = "fixed",
    stringsAsFactors = FALSE
  )
}

select_ar1_target <- function(scores, target = AR1_NULL_TARGET) {
  admissible <- scores[scores$null_f1 >= target, , drop = FALSE]
  if (!nrow(admissible)) {
    stop("AR(1) calibration grid does not attain null F1 target ", target)
  }
  admissible[which.min(admissible$candidate), , drop = FALSE]
}

select_ar1_closest <- function(scores, target) {
  scores[which.min(abs(scores$null_f1 - target)), , drop = FALSE]
}

attach_ar1_metadata <- function(
    scores, stage, target, selected, calibration_reps, validation_reps,
    calibration_seed, validation_seed) {
  scores$stage <- stage
  scores$target_f1 <- target
  scores$target_difference <- scores$null_f1 - target
  scores$selected <- (
    scores$method == "DeCAFS AR1" & selected$decafs_constant == scores$candidate
  ) | (
    scores$method == "PELT AR1 approximate" &
      selected$approximate_constant == scores$candidate
  ) | (
    scores$method == "SVP AR1Focus" &
      selected$svp_constant == scores$candidate
  )
  scores$calibration_reps <- calibration_reps
  scores$validation_reps <- validation_reps
  scores$calibration_seed <- calibration_seed
  scores$validation_seed <- validation_seed
  scores$target_nominal <- AR1_NULL_TARGET
  scores$decafs_version <- decafs_version()
  scores
}

ar1_method_labels <- function(calibration) {
  c(
    "PELT AR1 approximate",
    "PELT inflated",
    "DeCAFS AR1",
    paste0(
      "SVP AR1Focus + AR1 cost / right / c = ",
      formatC(calibration$selected$svp_constant, digits = 3L, format = "f")
    )
  )
}

run_ar1_calibration <- function(
    n = AR1_N,
    rho = AR1_RHO,
    calibration_reps = AR1_CALIBRATION_REPS,
    validation_reps = AR1_VALIDATION_REPS,
    calibration_seed = AR1_CALIBRATION_SEED,
    validation_seed = AR1_VALIDATION_SEED,
    decafs_grid = seq(1, 8, 0.25),
    approximate_grid = seq(0.5, 4, 0.05),
    svp_grid = seq(1, 8, 0.05),
    workers = 1L,
    write_outputs = TRUE) {
  require_decafs()
  calibration_null <- ar1_null_stream(
    n, rho, calibration_reps, calibration_seed
  )
    decafs_scores <- evaluate_ar1_candidates(
    decafs_grid, "DeCAFS AR1", calibration_null, n, rho, workers = workers
  )
  approximate_scores <- evaluate_ar1_candidates(
    approximate_grid, "PELT AR1 approximate", calibration_null, n, rho,
    workers = workers
  )
  svp_scores <- evaluate_ar1_candidates(
    svp_grid, "SVP AR1Focus", calibration_null, n, rho, workers = workers
  )
  selected_decafs <- select_ar1_target(decafs_scores)
  target_f1 <- selected_decafs$null_f1
  selected_approximate <- select_ar1_closest(approximate_scores, target_f1)
  selected_svp <- select_ar1_closest(svp_scores, target_f1)
  selected <- list(
    decafs_constant = selected_decafs$candidate,
    approximate_constant = selected_approximate$candidate,
    svp_constant = selected_svp$candidate,
    target_f1 = target_f1,
    target_f1_lower = selected_decafs$null_f1_lower,
    target_f1_upper = selected_decafs$null_f1_upper,
    n = n,
    rho = rho,
    tolerance = AR1_TOLERANCE,
    calibration_reps = calibration_reps,
    validation_reps = validation_reps,
    calibration_seed = calibration_seed,
    validation_seed = validation_seed,
    decafs_version = decafs_version()
  )

  validation_null <- ar1_null_stream(n, rho, validation_reps, validation_seed)
  validation <- dplyr::bind_rows(
    evaluate_ar1_candidates(
      selected_decafs$candidate, "DeCAFS AR1", validation_null, n, rho,
      workers = workers
    ),
    evaluate_ar1_candidates(
      selected_approximate$candidate, "PELT AR1 approximate",
      validation_null, n, rho, workers = workers
    ),
    evaluate_ar1_candidates(
      selected_svp$candidate, "SVP AR1Focus", validation_null, n, rho,
      workers = workers
    ),
    evaluate_fixed_ar1_method(
      "PELT inflated", validation_null, n, rho, validation_seed
    )
  )
  calibration_table <- dplyr::bind_rows(
    attach_ar1_metadata(
      decafs_scores, "calibration", target_f1, selected,
      calibration_reps, validation_reps, calibration_seed, validation_seed
    ),
    attach_ar1_metadata(
      approximate_scores, "calibration", target_f1, selected,
      calibration_reps, validation_reps, calibration_seed, validation_seed
    ),
    attach_ar1_metadata(
      svp_scores, "calibration", target_f1, selected,
      calibration_reps, validation_reps, calibration_seed, validation_seed
    ),
    attach_ar1_metadata(
      validation, "validation", target_f1, selected,
      calibration_reps, validation_reps, calibration_seed, validation_seed
    ),
    evaluate_fixed_ar1_method(
      "PELT inflated", calibration_null, n, rho, calibration_seed
    )
  )
  calibration_table$target_nominal <- AR1_NULL_TARGET
  calibration_table$calibration_seed <- calibration_seed
  calibration_table$validation_seed <- validation_seed
  result <- list(
    selected = selected,
    target_nominal = AR1_NULL_TARGET,
    decafs_version = decafs_version(),
    candidate_tables = list(
      decafs = decafs_scores,
      approximate = approximate_scores,
      svp = svp_scores
    ),
    table = calibration_table,
    validation = validation,
    fixed_pelt = list(
      approximate = list(constant = 2.2),
      inflated = list(constant = 3, inflation = (1 + rho) / (1 - rho))
    )
  )
  if (write_outputs) {
    dir.create(AR1_ROOT, recursive = TRUE, showWarnings = FALSE)
    saveRDS(result, file.path(AR1_ROOT, "ar1_calibration.rds"))
    utils::write.csv(
      calibration_table,
      file.path(AR1_ROOT, "true_true_calibration.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      dplyr::bind_rows(decafs_scores, approximate_scores, svp_scores),
      file.path(AR1_ROOT, "ar1_null_calibration.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      data.frame(
        method = c("DeCAFS AR1", "PELT AR1 approximate", "SVP AR1Focus"),
        parameter = c(
          selected$decafs_constant,
          selected$approximate_constant,
          selected$svp_constant
        ),
        target_f1 = target_f1,
        decafs_version = decafs_version(),
        stringsAsFactors = FALSE
      ),
      file.path(AR1_ROOT, "selected_parameters.csv"),
      row.names = FALSE
    )
  }
  result
}

if (identical(tolower(Sys.getenv("SVP_RUN_AR1_CALIBRATION")), "true")) {
  print(run_ar1_calibration())
}
