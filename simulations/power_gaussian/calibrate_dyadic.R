## Jointly tune dyadic SVP thresholds under a 1% no-change FPR constraint.
library(svpChange)
source(file.path("simulations", "power_common.R"))
source(file.path("simulations", "gaussian_common.R"))

GAUSSIAN_DYADIC_TUNING_FILE <- file.path(
  "simulations", "power_gaussian", "gaussian_dyadic_calibration.csv"
)
GAUSSIAN_DYADIC_SEARCH_FILE <- file.path(
  "simulations", "power_gaussian", "gaussian_dyadic_search.csv"
)

dyadic_wilson_interval <- function(successes, trials, level = 0.95) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  p <- successes / trials
  denominator <- 1 + z^2 / trials
  center <- (p + z^2 / (2 * trials)) / denominator
  half_width <- z * sqrt(p * (1 - p) / trials + z^2 / (4 * trials^2)) /
    denominator
  c(lower = max(0, center - half_width),
    upper = min(1, center + half_width))
}

dyadic_n_bands <- function(n) max(1L, floor(log2(n / 2)) + 1L)

dyadic_band_limits <- function(n) {
  bands <- dyadic_n_bands(n)
  data.frame(
    band = 0:(bands - 1L),
    min_length = vapply(0:(bands - 1L), function(j) {
      as.integer(floor(n / 2^(j + 1L)) + 1L)
    }, integer(1)),
    max_length = vapply(0:(bands - 1L), function(j) {
      as.integer(floor(n / 2^j))
    }, integer(1))
  )
}

dyadic_null_false_positive <- function(i, n, seed, gamma) {
  set_power_seed(seed + i)
  y <- stats::rnorm(n)
  length(SVP_both_dyadic(y, gamma)$changepoints) > 1L
}

dyadic_updown_f1 <- function(i, n, jump_sizes, reps, seed, gamma,
                             tolerance) {
  jump_index <- (i - 1L) %/% reps + 1L
  rep_index <- (i - 1L) %% reps + 1L
  pattern_index <- match("updown", POWER_PATTERNS)
  design_index <- pattern_index + (jump_index - 1L) *
    length(POWER_PATTERNS) + (rep_index - 1L) * length(POWER_PATTERNS) *
    length(jump_sizes)
  set_power_seed(seed + design_index)
  mu <- generate_power_signal(n, "updown", jump_sizes[jump_index])
  y <- mu + stats::rnorm(n)
  truth <- normalise_boundaries(which(diff(mu) != 0), n)
  estimate <- normalise_boundaries(
    SVP_both_dyadic(y, gamma)$changepoints, n
  )
  paper_metrics(truth, estimate, tolerance, terminal = n)[["F1"]]
}

score_dyadic_vector <- function(
    gamma, n, jump_sizes, tuning_reps, calibration_reps,
    calibration_seed, tuning_seed, workers, tolerance,
    target = 0.01) {
  false_positive <- power_lapply(seq_len(calibration_reps), function(i) {
    dyadic_null_false_positive(i, n, calibration_seed, gamma)
  }, workers = workers)
  false_detection <- mean(unlist(false_positive, use.names = FALSE))
  mean_f1 <- NA_real_
  if (false_detection <= target) {
    n_tuning <- length(jump_sizes) * tuning_reps
    f1 <- power_lapply(seq_len(n_tuning), function(i) {
      dyadic_updown_f1(i, n, jump_sizes, tuning_reps, tuning_seed,
                       gamma, tolerance)
    }, workers = workers)
    mean_f1 <- mean(unlist(f1, use.names = FALSE))
  }
  data.frame(
    null_calibration_fpr = false_detection,
    mean_updown_f1 = mean_f1,
    stringsAsFactors = FALSE
  )
}

select_dyadic_candidate <- function(rows, gamma_vectors, target = 0.01) {
  feasible <- which(is.finite(rows$mean_updown_f1) &
                      rows$null_calibration_fpr <= target)
  if (!length(feasible)) {
    best_fpr <- min(rows$null_calibration_fpr)
    closest <- which(rows$null_calibration_fpr == best_fpr)
    vectors <- gamma_vectors[closest]
    lex <- do.call(rbind, vectors)
    order_args <- lapply(seq_len(ncol(lex)), function(j) -lex[, j])
    return(closest[do.call(order, order_args)[1L]])
  }
  values <- rows[feasible, , drop = FALSE]
  vectors <- gamma_vectors[feasible]
  lex <- do.call(rbind, vectors)
  names <- paste0("g", seq_len(ncol(lex)))
  colnames(lex) <- names
  order_args <- c(
    list(-values$mean_updown_f1,
         abs(values$null_calibration_fpr - target)),
    lapply(seq_len(ncol(lex)), function(j) -lex[, j])
  )
  feasible[do.call(order, order_args)[1L]]
}

run_gaussian_dyadic_calibration <- function(
    n = 1000L, jump_sizes = seq(0.1, 2, 0.1), tuning_reps = 100L,
    calibration_reps = 20000L, validation_reps = 10000L,
    search_tuning_reps = 10L, search_calibration_reps = 500L,
    calibration_seed = 2100000L, validation_seed = 3000000L,
    search_seed = 2000000L,
    tuning_seed = 1000000L, workers = power_default_workers(),
    max_sweeps = 5L, n_finalists = 5L, target = 0.01) {
  if (n < 2L || n != as.integer(n)) stop("n must be an integer of at least 2")
  if (tuning_reps < 1L || calibration_reps < 1L || validation_reps < 1L ||
      search_tuning_reps < 1L || search_calibration_reps < 1L)
    stop("replication counts must be positive")
  if (max_sweeps < 1L) stop("max_sweeps must be positive")
  if (n_finalists < 1L) stop("n_finalists must be positive")

  bands <- dyadic_n_bands(n)
  reference_gamma <- read_gaussian_calibration()[["both"]] * log(n)
  tolerance <- round(n * 0.0025)
  jump_sizes <- as.numeric(jump_sizes)
  current <- rep(reference_gamma, bands)
  history <- list()

  evaluate_vectors <- function(vectors, phase, sweep, band) {
    message("Scoring ", length(vectors), " ", phase,
            " candidates at band ", band, ", sweep ", sweep)
    scores <- power_lapply(seq_along(vectors), function(i) {
      score_dyadic_vector(
        vectors[[i]], n, jump_sizes, search_tuning_reps,
        search_calibration_reps, search_seed, tuning_seed, 1L,
        tolerance, target
      )
    }, workers = min(workers, length(vectors)))
    rows <- do.call(rbind, scores)
    rows$phase <- phase
    rows$sweep <- sweep
    rows$band <- band
    rows$gamma <- vapply(vectors, function(value) {
      paste(format(value, digits = 17, scientific = FALSE, trim = TRUE),
            collapse = ";")
    }, character(1))
    history[[length(history) + 1L]] <<- rows
    rows
  }

  coordinate_search <- function(multipliers, phase, local = FALSE) {
    for (sweep in seq_len(max_sweeps)) {
      changed <- FALSE
      for (band in seq_len(bands)) {
        band_multipliers <- if (local) {
          center <- current[band] / reference_gamma
          seq(max(0.5, center - 0.25), min(3.5, center + 0.25),
              by = 0.05)
        } else {
          multipliers
        }
        values <- sort(unique(c(reference_gamma * band_multipliers,
                                current[band])))
        vectors <- lapply(values, function(value) {
          candidate <- current
          candidate[band] <- value
          candidate
        })
        rows <- evaluate_vectors(vectors, phase, sweep, band - 1L)
        selected <- select_dyadic_candidate(rows, vectors, target)
        if (is.na(selected)) next
        if (!identical(as.numeric(vectors[[selected]]),
                       as.numeric(current))) {
          current <<- vectors[[selected]]
          changed <- TRUE
        }
      }
      if (!changed) break
    }
    current
  }

  coordinate_search(seq(0.5, 3.5, by = 0.25), "coarse")
  coordinate_search(NULL, "refine", local = TRUE)

  screen_history <- do.call(rbind, history)
  screened <- which(is.finite(screen_history$mean_updown_f1) &
                      screen_history$null_calibration_fpr <= target)
  if (length(screened)) {
    screen_order <- screened[order(
      -screen_history$mean_updown_f1[screened],
      abs(screen_history$null_calibration_fpr[screened] - target)
    )]
    parse_vector <- function(value) {
      as.numeric(strsplit(value, ";", fixed = TRUE)[[1L]])
    }
    finalists <- lapply(head(screen_order, n_finalists), function(i) {
      parse_vector(screen_history$gamma[i])
    })
  } else {
    finalists <- list()
  }
  finalists <- c(finalists, list(current, rep(reference_gamma, bands),
                                 rep(1e6, bands)))
  finalist_keys <- vapply(finalists, function(value) {
    paste(format(value, digits = 17, scientific = FALSE, trim = TRUE),
          collapse = ";")
  }, character(1))
  finalists <- finalists[!duplicated(finalist_keys)]

  message("Re-scoring ", length(finalists),
          " finalists on the full calibration and tuning samples")
  finalist_scores <- lapply(finalists, function(value) {
    score_dyadic_vector(
      value, n, jump_sizes, tuning_reps, calibration_reps,
      calibration_seed, tuning_seed, workers, tolerance, target
    )
  })
  finalist_scores <- do.call(rbind, finalist_scores)
  selected_index <- select_dyadic_candidate(finalist_scores, finalists, target)
  if (is.na(selected_index) ||
      !is.finite(finalist_scores$mean_updown_f1[selected_index])) {
    stop("No dyadic finalist meets the full-sample null FPR target")
  }
  current <- finalists[[selected_index]]
  final_score <- finalist_scores[selected_index, , drop = FALSE]

  validation <- power_lapply(seq_len(validation_reps), function(i) {
    dyadic_null_false_positive(i, n, validation_seed, current)
  }, workers = workers)
  validation_count <- sum(unlist(validation, use.names = FALSE))
  validation_fpr <- validation_count / validation_reps
  validation_interval <- dyadic_wilson_interval(validation_count,
                                                  validation_reps)

  limits <- dyadic_band_limits(n)
  selected <- data.frame(
    limits,
    gamma = current,
    selected = TRUE,
    n = n,
    null_calibration_fpr = final_score$null_calibration_fpr,
    null_validation_fpr = validation_fpr,
    validation_lower = validation_interval[["lower"]],
    validation_upper = validation_interval[["upper"]],
    target_fpr = target,
    mean_updown_f1 = final_score$mean_updown_f1,
    tuning_reps_per_jump = tuning_reps,
    tuning_jump_count = length(jump_sizes),
    search_tuning_reps_per_jump = search_tuning_reps,
    search_calibration_reps = search_calibration_reps,
    n_finalists = length(finalists),
    calibration_reps = calibration_reps,
    validation_reps = validation_reps,
    calibration_seed = calibration_seed,
    validation_seed = validation_seed,
    tuning_seed = tuning_seed,
    search = "screened coordinate coarse 0.25, local refine 0.05; multipliers 0.5-3.5",
    stringsAsFactors = FALSE
  )
  history <- do.call(rbind, history)
  utils::write.csv(selected, GAUSSIAN_DYADIC_TUNING_FILE, row.names = FALSE)
  utils::write.csv(history, GAUSSIAN_DYADIC_SEARCH_FILE, row.names = FALSE)
  raw <- list(selected = selected, history = history,
              finalist_scores = finalist_scores,
              finalists = finalists,
              gamma = current, final_score = final_score,
              validation_false_positives = validation_count,
              validation_interval = validation_interval,
              n = n, jump_sizes = jump_sizes,
              tuning_reps = tuning_reps,
              calibration_reps = calibration_reps,
              validation_reps = validation_reps,
              search_tuning_reps = search_tuning_reps,
              search_calibration_reps = search_calibration_reps,
              calibration_seed = calibration_seed,
              validation_seed = validation_seed, search_seed = search_seed,
              tuning_seed = tuning_seed,
              search = list(coarse = seq(0.5, 3.5, by = 0.25),
                            refine_step = 0.05, refine_radius = 0.25,
                            max_sweeps = max_sweeps))
  saveRDS(raw, sub("\\.csv$", ".rds", GAUSSIAN_DYADIC_TUNING_FILE))
  message("Selected dyadic gamma: ", paste(format(current, digits = 6),
                                             collapse = ", "),
          "; calibration FPR = ", format(final_score$null_calibration_fpr,
                                            digits = 4),
          "; validation FPR = ", format(validation_fpr, digits = 4),
          "; mean updown F1 = ", format(final_score$mean_updown_f1,
                                           digits = 5))
  selected
}

if (identical(tolower(Sys.getenv("SVP_CALIBRATE_GAUSSIAN_DYADIC")), "true"))
  print(run_gaussian_dyadic_calibration(workers = 8L))
