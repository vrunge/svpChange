## Benchmark stepR::stepFit() against svpChange::svp_smuce_cpp().
## Run from the package root after installing stepR.

if (!requireNamespace("stepR", quietly = TRUE))
  stop("Install stepR first with install.packages('stepR')")
if (!requireNamespace("svpChange", quietly = TRUE))
  stop("Install svpChange first")

benchmark_smuce_time <- function(
    n_values = c(50L, 100L, 200L, 400L, 800L),
    repetitions = 3L,
    sigma = 1,
    q = 0.8,
    seed = NULL) {
  stopifnot(sigma > 0, q >= 0)
  if (!is.null(seed)) set.seed(seed)
  out <- vector("list", length(n_values) * repetitions)
  z <- 1L
  for (n in n_values) for (replicate in seq_len(repetitions)) {
    y <- as.double(rnorm(n, sd = sigma))
    step_call <- function() stepR::stepFit(
      y = y, q = q, family = "gauss", sd = sigma,
      intervalSystem = "all", lengths = seq_len(n), penalty = "sqrt",
      confband = FALSE, jumpint = FALSE
    )
    svp_call <- function() svpChange::svp_smuce_cpp(y, q, sigma^2)

    ## Warm up and verify equality before recording timings.
    step_fit <- step_call()
    svp_fit <- svp_call()
    stopifnot(identical(as.integer(step_fit$rightIndex), svp_fit))

    smuce <- system.time(step_call())
    cpp <- system.time(svp_call())
    out[[z]] <- data.frame(
      n = n, replicate = replicate,
      smuce_seconds = unname(smuce[["elapsed"]]),
      svp_smuce_cpp_seconds = unname(cpp[["elapsed"]])
    )
    z <- z + 1L
  }
  do.call(rbind, out)
}

## Example: timings <- benchmark_smuce_time(seed = NULL); print(timings)
