## Audit the refreshed Gaussian timing results against the two preserved
## historical result sets. Absolute ratios are descriptive because the older
## files do not record a complete, matching machine/timer manifest.

library(dplyr)

TIME_ROOT <- file.path("simulations", "time_gaussian")
LEGACY_ROOT <- file.path("simulations", "other_simus", "time_gaussian_legacy")

timing_slope_summary <- function(data, dataset) {
  data |>
    filter(experiment == "vary_n", is.finite(time), time > 0) |>
    group_by(method, n) |>
    summarise(median_time = median(time), .groups = "drop") |>
    group_by(method) |>
    summarise(
      dataset = dataset,
      exponent = unname(coef(lm(log(median_time) ~ log(n)))[2L]),
      n_min = min(n), n_max = max(n),
      median_time_at_n_max = median_time[which.max(n)],
      .groups = "drop"
    ) |>
    select(dataset, everything())
}

read_wide_legacy <- function(path) {
  data <- read.csv(path, stringsAsFactors = FALSE)
  bind_rows(
    transmute(data, experiment = "vary_n", n, method = "PELT",
              time = PELT),
    transmute(data, experiment = "vary_n", n, method = "SVP historical",
              time = SVP)
  )
}

compare_gaussian_time <- function(root = TIME_ROOT) {
  current <- read.csv(
    file.path(root, "gaussian_time_results.csv"), stringsAsFactors = FALSE
  )
  legacy <- read_wide_legacy(file.path(
    LEGACY_ROOT, "complexity_exp1_nochange.csv"
  ))
  paper_path <- file.path(
    "..", "svpChange_old", "simulations_paper", "time_gaussian",
    "gaussian_time_results.csv"
  )
  summaries <- list(
    timing_slope_summary(legacy, "historical Experiment 1 CSV"),
    timing_slope_summary(current, "refreshed controlled benchmark")
  )
  if (file.exists(paper_path)) {
    paper <- read.csv(paper_path, stringsAsFactors = FALSE)
    old_svp_names <- setdiff(unique(paper$method), "PELT")
    paper$method[paper$method != "PELT"] <- paste(
      "SVP historical curve",
      LETTERS[match(paper$method[paper$method != "PELT"], old_svp_names)]
    )
    summaries <- append(
      summaries,
      list(timing_slope_summary(paper, "checked-in paper timing CSV")),
      after = 1L
    )
  }
  comparison <- bind_rows(summaries)
  write.csv(
    comparison, file.path(root, "gaussian_time_historical_comparison.csv"),
    row.names = FALSE
  )

  paper_script <- file.path(
    "..", "svpChange_old", "simulations_paper", "time_gaussian", "run_time.R"
  )
  notes <- c(
    "# Gaussian timing audit",
    "",
    "- The historical Experiment 1 calls `SVP()` without an explicit mode;",
    "  its meaning depends on the package version installed when it was run.",
    "- The checked-in paper timing script currently assigns the same `both`,",
    "  `c = 1.5` call to two differently named curves. Its saved CSV does not",
    "  show identical timings for those curves, so script and result provenance",
    "  are not internally reproducible as checked in.",
    "- The refreshed script specifies every mode and constant, warms each method,",
    "  randomizes method order, shares data, and batches sub-resolution calls.",
    "- Absolute old/new time ratios are not treated as hardware-independent.",
    paste0("- Audited paper script: `", paper_script, "`.")
  )
  writeLines(notes, file.path(root, "gaussian_time_historical_comparison.md"))
  invisible(comparison)
}

if (identical(tolower(Sys.getenv("SVP_COMPARE_GAUSSIAN_TIME")), "true")) {
  compare_gaussian_time()
}
