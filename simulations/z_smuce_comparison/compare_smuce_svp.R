## Compare stepR::stepFit() with the multiscale Gaussian SVP routines.

if (!requireNamespace("stepR", quietly = TRUE)) {
  stop("Install stepR first with install.packages('stepR')")
}
if (!requireNamespace("svpChange", quietly = TRUE)) {
  stop("Install svpChange first")
}

library(svpChange)

compare_smuce_svp <- function(
    n = 100,
    sigma = 0.5,
    q = 0.8,
    seed = NULL,
    n_segments = 10) {

  stopifnot(n >= 10, sigma > 0, q >= 0)
  if (!is.null(seed)) set.seed(seed)

  ## Random heterogeneous piecewise-constant Gaussian signal.
  cuts <- sort(sample(2:(n - 1), n_segments - 1L))
  ends <- c(cuts, n)
  lengths <- diff(c(0, ends))
  means <- cumsum(stats::rnorm(n_segments, mean = 0, sd = 1.5))
  signal <- rep(means, lengths)
  y <- as.double(signal + stats::rnorm(n, sd = sigma))

  smuce_fit <- stepR::stepFit(
    y = y,
    q = q,
    family = "gauss",
    sd = sigma,
    intervalSystem = "all",
    lengths = seq_len(n),
    penalty = "sqrt",
    confband = FALSE,
    jumpint = FALSE
  )

  constrained <- svpChange::svp_smuce(y, q, sigma^2)
  constrained_cpp <- svpChange::svp_smuce_cpp(y, q, sigma^2)

  list(
    smuce = as.integer(smuce_fit$rightIndex),
    svp_smuce = constrained,
    svp_smuce_cpp = constrained_cpp
  )
}

#' Run repeated random SMUCE/SVP comparisons.
#' @return A data.frame with one row per simulation and an equality indicator.
run_smuce_replicates <- function(B = 20, n = 120, ...) {
  out <- vector("list", B)
  for (b in seq_len(B)) {
    z <- compare_smuce_svp(n = n, seed = NULL, ...)
    out[[b]] <- data.frame(
      replicate = b,
      equal = identical(z$smuce, z$svp_smuce),
      equal_cpp = identical(z$svp_smuce, z$svp_smuce_cpp),
      smuce = paste(z$smuce, collapse = ","),
      svp_smuce = paste(z$svp_smuce, collapse = ","),
      svp_smuce_cpp = paste(z$svp_smuce_cpp, collapse = ","),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

## Run explicitly with:
## source("simulations/smuce_comparison/compare_smuce_svp.R")
## answer <- compare_smuce_svp(seed = 1)
## answer
