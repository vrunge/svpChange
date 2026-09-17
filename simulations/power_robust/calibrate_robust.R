## Independent null calibration for the Student-t(2) robust power study.
## Run from the svpChange package root.

library(svpChange)
library(changepoint)
library(robseg)
library(dplyr)

source(file.path("simulations", "power_common.R"))

ROBUST_CALIBRATION_N <- 1000L
ROBUST_CALIBRATION_REPS <- 20000L
ROBUST_VALIDATION_REPS <- 10000L
ROBUST_CALIBRATION_SEED <- 910000L
ROBUST_VALIDATION_SEED <- 920000L
ROBUST_RFPOP_TARGET <- 0.99
ROBUST_CALIBRATION_ROOT <- file.path("simulations", "power_robust")
## robseg names its biweight loss "Outlier".
RFPOP_LOSS <- "Outlier"

student_null <- function(n, seed) {
  set_power_seed(seed)
  ts_generator(
    chpts = n, parameters = 0, type = "student", df = 2, scale = 1
  )
}

student_null_stream <- function(n, reps, seed) {
  lapply(seq_len(reps), function(i) student_null(n, seed + i))
}

wilson_interval <- function(successes, trials, level = 0.95) {
  if (trials < 1L) stop("trials must be positive")
  z <- stats::qnorm(1 - (1 - level) / 2)
  p <- successes / trials
  denominator <- 1 + z^2 / trials
  centre <- (p + z^2 / (2 * trials)) / denominator
  half_width <- z * sqrt(p * (1 - p) / trials + z^2 / (4 * trials^2)) /
    denominator
  c(lower = max(0, centre - half_width),
    upper = min(1, centre + half_width))
}

fit_rfpop <- function(y, penalty_constant) {
  scale <- stats::mad(y, constant = 1.4826)
  if (!is.finite(scale) || scale <= 0) {
    stop("RFPOP requires a positive finite MAD scale")
  }
  fit <- robseg::Rob_seg.std(
    y / scale, loss = RFPOP_LOSS,
    lambda = penalty_constant * log(length(y)), lthreshold = 3
  )
  list(
    boundaries = normalise_boundaries(head(fit$t.est, -1L), length(y)),
    fitted = as.numeric(fit$smt) * scale
  )
}

mood_threshold <- function(n, oracle_segments, alpha) {
  if (length(n) != 1L || n < 2L || n %% 1 != 0) {
    stop("n must be one integer greater than one")
  }
  if (length(oracle_segments) != 1L || oracle_segments < 1L ||
      oracle_segments > n || oracle_segments %% 1 != 0) {
    stop("oracle_segments must be one integer in [1, n]")
  }
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be one finite value in (0, 1)")
  }
  segment_length <- n / oracle_segments
  if (segment_length <= 1) stop("oracle segments leave no valid split")
  alpha_s <- 1 - (1 - alpha)^(1 / (segment_length - 1))
  stats::qchisq(1 - alpha_s, df = 1)
}

wilcoxon_threshold <- function(n, oracle_segments, constant) {
  if (length(constant) != 1L || !is.finite(constant) || constant <= 0) {
    stop("constant must be one finite positive value")
  }
  constant * sqrt((n / oracle_segments)^3 / 12)
}

null_f1 <- function(boundaries, n, tolerance = round(n * 0.0025)) {
  paper_metrics(c(n), boundaries, tolerance, terminal = n)[["F1"]]
}

evaluate_null_candidates <- function(
    candidates, method, n, reps, seed, oracle_segments = 1L,
    workers = 1L, null_series = NULL) {
  if (is.null(null_series)) {
    null_series <- student_null_stream(n, reps, seed)
  }
  if (length(null_series) != reps) stop("null_series has the wrong length")
  score_one <- function(candidate) {
    successes <- logical(reps)
    for (i in seq_len(reps)) {
      y <- null_series[[i]]
      boundaries <- switch(
        method,
        RFPOP = fit_rfpop(y, candidate)$boundaries,
        Wilcoxon = normalise_boundaries(
          SVP(
            y, wilcoxon_threshold(n, oracle_segments, candidate),
            "WilcoxonCost", subtests = "right"
          )$changepoints, n
        ),
        MedianMood = normalise_boundaries(
          SVP(
            y, mood_threshold(n, oracle_segments, candidate),
            "MedianMoodCost", subtests = "right"
          )$changepoints, n
        ),
        stop("Unknown calibration method: ", method)
      )
      successes[i] <- null_f1(boundaries, n) == 1
    }
    interval <- wilson_interval(sum(successes), reps)
    data.frame(
      method = method,
      candidate = candidate,
      null_f1 = mean(successes),
      null_f1_lower = interval[["lower"]],
      null_f1_upper = interval[["upper"]],
      recovered_one_segment = sum(successes),
      reps = reps,
      seed = seed,
      n = n,
      oracle_segments = oracle_segments,
      stringsAsFactors = FALSE
    )
  }
  power_lapply(candidates, score_one, workers = workers) |>
    dplyr::bind_rows()
}

validate_candidate <- function(
    method, candidate, n, reps, seed, oracle_segments = 1L,
    workers = 1L, null_series = NULL) {
  evaluate_null_candidates(
    candidate, method, n, reps, seed, oracle_segments, workers, null_series
  )
}

pelt_null_benchmark <- function(n, reps, seed, workers = 1L) {
  score_one <- function(i) {
    y <- student_null(n, seed + i)
    fit <- changepoint::cpt.mean(
      y, method = "PELT", penalty = "Manual", pen.value = 2 * log(n)
    )
    null_f1(normalise_boundaries(changepoint::cpts(fit), n), n)
  }
  scores <- unlist(power_lapply(seq_len(reps), score_one, workers = workers))
  interval <- wilson_interval(sum(scores == 1), reps)
  data.frame(
    student_pelt_null_f1 = mean(scores),
    student_pelt_null_f1_lower = interval[["lower"]],
    student_pelt_null_f1_upper = interval[["upper"]],
    student_pelt_reps = reps,
    student_pelt_seed = seed,
    stringsAsFactors = FALSE
  )
}

gaussian_pelt_benchmark <- function() {
  results_path <- file.path("simulations", "power_gaussian", "results.csv")
  if (file.exists(results_path)) {
    results <- utils::read.csv(results_path, stringsAsFactors = FALSE)
    scores <- results$F1[
      results$algorithm == "PELT" & results$pattern == "none"
    ]
    interval <- wilson_interval(sum(scores == 1), length(scores))
    return(data.frame(
      gaussian_pelt_null_f1 = mean(scores),
      gaussian_pelt_null_f1_lower = interval[["lower"]],
      gaussian_pelt_null_f1_upper = interval[["upper"]],
      gaussian_pelt_reps = length(scores),
      stringsAsFactors = FALSE
    ))
  }
  path <- file.path("simulations", "power_gaussian",
                    "gaussian_mode_calibration.csv")
  if (!file.exists(path)) {
    return(data.frame(
      gaussian_pelt_null_f1 = NA_real_,
      gaussian_pelt_null_f1_lower = NA_real_,
      gaussian_pelt_null_f1_upper = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  calibration <- utils::read.csv(path, stringsAsFactors = FALSE)
  row <- calibration[calibration$mode == "right", ][1L, ]
  data.frame(
    gaussian_pelt_null_f1 = 1 - row$pelt_fpr_validation,
    gaussian_pelt_null_f1_lower = 1 - row$pelt_fpr_upper,
    gaussian_pelt_null_f1_upper = 1 - row$pelt_fpr_lower,
    stringsAsFactors = FALSE
  )
}

attach_calibration_metadata <- function(
    scores, stage, target, target_interval, validation_reps,
    validation_seed, selected, gaussian_benchmark, student_benchmark) {
  scores$stage <- stage
  scores$rfpop_target_f1 <- target
  scores$rfpop_target_lower <- target_interval[["lower"]]
  scores$rfpop_target_upper <- target_interval[["upper"]]
  scores$target_difference <- scores$null_f1 - target
  scores$selected <-
    scores$method == "RFPOP" & selected$rfpop_constant == scores$candidate |
    scores$method == "Wilcoxon" & selected$wilcoxon_constant == scores$candidate |
    scores$method == "MedianMood" & selected$mood_alpha == scores$candidate
  scores$validation_reps <- validation_reps
  scores$validation_seed <- validation_seed
  scores$gaussian_pelt_null_f1 <- gaussian_benchmark$gaussian_pelt_null_f1
  scores$gaussian_pelt_null_f1_lower <-
    gaussian_benchmark$gaussian_pelt_null_f1_lower
  scores$gaussian_pelt_null_f1_upper <-
    gaussian_benchmark$gaussian_pelt_null_f1_upper
  scores$student_pelt_null_f1 <- student_benchmark$student_pelt_null_f1
  scores$student_pelt_null_f1_lower <-
    student_benchmark$student_pelt_null_f1_lower
  scores$student_pelt_null_f1_upper <-
    student_benchmark$student_pelt_null_f1_upper
  scores
}

run_robust_calibration <- function(
    n = ROBUST_CALIBRATION_N,
    calibration_reps = ROBUST_CALIBRATION_REPS,
    validation_reps = ROBUST_VALIDATION_REPS,
    calibration_seed = ROBUST_CALIBRATION_SEED,
    validation_seed = ROBUST_VALIDATION_SEED,
    rfpop_target = ROBUST_RFPOP_TARGET,
    rfpop_grid = seq(2, 4, 0.5),
    wilcoxon_grid = seq(1.25, 2.25, 0.25),
    mood_grid = c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.02),
    workers = 1L, write_outputs = TRUE) {
  calibration_null <- student_null_stream(n, calibration_reps, calibration_seed)
  rfpop_calibration <- evaluate_null_candidates(
    rfpop_grid, "RFPOP", n, calibration_reps, calibration_seed,
    oracle_segments = 1L, workers = workers, null_series = calibration_null
  )
  admissible <- rfpop_calibration[rfpop_calibration$null_f1 >= rfpop_target, ]
  if (!nrow(admissible)) {
    stop("RFPOP calibration grid does not attain the nominal null F1 target")
  }
  selected_rfpop <- admissible[which.min(admissible$candidate), ]
  selected_rfpop_constant <- selected_rfpop$candidate
  rfpop_target_interval <- c(
    lower = selected_rfpop$null_f1_lower,
    upper = selected_rfpop$null_f1_upper
  )

  wilcoxon_calibration <- evaluate_null_candidates(
    wilcoxon_grid, "Wilcoxon", n, calibration_reps, calibration_seed,
    oracle_segments = 1L, workers = workers, null_series = calibration_null
  )
  mood_calibration <- evaluate_null_candidates(
    mood_grid, "MedianMood", n, calibration_reps, calibration_seed,
    oracle_segments = 1L, workers = workers, null_series = calibration_null
  )
  select_closest <- function(scores) {
    scores[which.min(abs(scores$null_f1 - selected_rfpop$null_f1)), ]
  }
  selected_wilcoxon <- select_closest(wilcoxon_calibration)
  selected_mood <- select_closest(mood_calibration)
  selected <- list(
    rfpop_constant = selected_rfpop_constant,
    wilcoxon_constant = selected_wilcoxon$candidate,
    mood_alpha = selected_mood$candidate,
    rfpop_target_f1 = selected_rfpop$null_f1,
    rfpop_target_f1_lower = selected_rfpop$null_f1_lower,
    rfpop_target_f1_upper = selected_rfpop$null_f1_upper,
    n = n, oracle_segments_null = 1L,
    calibration_reps = calibration_reps,
    validation_reps = validation_reps,
    calibration_seed = calibration_seed,
    validation_seed = validation_seed
  )

  validation_null <- student_null_stream(n, validation_reps, validation_seed)
  validation <- dplyr::bind_rows(
    validate_candidate(
      "RFPOP", selected_rfpop_constant, n, validation_reps,
      validation_seed, oracle_segments = 1L, workers = workers,
      null_series = validation_null
    ),
    validate_candidate(
      "Wilcoxon", selected_wilcoxon$candidate, n, validation_reps,
      validation_seed, oracle_segments = 1L, workers = workers,
      null_series = validation_null
    ),
    validate_candidate(
      "MedianMood", selected_mood$candidate, n, validation_reps,
      validation_seed, oracle_segments = 1L, workers = workers,
      null_series = validation_null
    )
  )
  student_scores <- vapply(validation_null, function(y) {
    fit <- changepoint::cpt.mean(
      y, method = "PELT", penalty = "Manual", pen.value = 2 * log(n)
    )
    null_f1(normalise_boundaries(changepoint::cpts(fit), n), n)
  }, numeric(1))
  student_interval <- wilson_interval(sum(student_scores == 1), validation_reps)
  student_benchmark <- data.frame(
    student_pelt_null_f1 = mean(student_scores),
    student_pelt_null_f1_lower = student_interval[["lower"]],
    student_pelt_null_f1_upper = student_interval[["upper"]],
    student_pelt_reps = validation_reps,
    student_pelt_seed = validation_seed,
    stringsAsFactors = FALSE
  )
  gaussian_benchmark <- gaussian_pelt_benchmark()
  all_scores <- dplyr::bind_rows(
    attach_calibration_metadata(
      rfpop_calibration, "calibration", selected_rfpop$null_f1,
      rfpop_target_interval, validation_reps, validation_seed, selected,
      gaussian_benchmark, student_benchmark
    ),
    attach_calibration_metadata(
      wilcoxon_calibration, "calibration", selected_rfpop$null_f1,
      rfpop_target_interval, validation_reps, validation_seed, selected,
      gaussian_benchmark, student_benchmark
    ),
    attach_calibration_metadata(
      mood_calibration, "calibration", selected_rfpop$null_f1,
      rfpop_target_interval, validation_reps, validation_seed, selected,
      gaussian_benchmark, student_benchmark
    ),
    attach_calibration_metadata(
      validation, "validation", selected_rfpop$null_f1,
      rfpop_target_interval, validation_reps, validation_seed, selected,
      gaussian_benchmark, student_benchmark
    )
  )
  all_scores$rfpop_target_nominal <- rfpop_target
  all_scores$calibration_seed <- calibration_seed
  result <- list(
    selected = selected,
    target_nominal = rfpop_target,
    table = all_scores,
    gaussian_pelt_benchmark = gaussian_benchmark,
    student_pelt_benchmark = student_benchmark,
    validation = validation
  )
  if (write_outputs) {
    dir.create(ROBUST_CALIBRATION_ROOT, recursive = TRUE,
               showWarnings = FALSE)
    saveRDS(result, file.path(ROBUST_CALIBRATION_ROOT,
                              "robust_calibration.rds"))
    utils::write.csv(
      all_scores,
      file.path(ROBUST_CALIBRATION_ROOT, "true_true_calibration.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      rfpop_calibration,
      file.path(ROBUST_CALIBRATION_ROOT, "rfpop_null_calibration.csv"),
      row.names = FALSE
    )
    selected_table <- data.frame(
      method = c("RFPOP", "Wilcoxon", "MedianMood"),
      parameter = c(selected$rfpop_constant, selected$wilcoxon_constant,
                    selected$mood_alpha),
      stringsAsFactors = FALSE
    )
    utils::write.csv(
      selected_table,
      file.path(ROBUST_CALIBRATION_ROOT, "selected_parameters.csv"),
      row.names = FALSE
    )
  }
  result
}
