.smuce_test_partitions <- function(n) {
  lapply(seq_len(2^(n - 1L)) - 1L, function(mask) {
    cuts <- which(as.logical(intToBits(mask)[seq_len(n - 1L)]))
    as.integer(c(cuts, n))
  })
}

.smuce_test_partition_cost <- function(y, ends, gamma, sigma2, n) {
  starts <- c(1L, head(ends, -1L) + 1L)
  sum(vapply(seq_along(ends), function(i) {
    smuce_cost(y[starts[i]:ends[i]], gamma, sigma2, n)
  }, numeric(1)))
}

test_that("R SMUCE dynamic programming is globally optimal", {
  set.seed(202)
  y <- c(rnorm(4), 3 + rnorm(4))
  gamma <- 1.5
  sigma2 <- 1
  n <- length(y)

  parts <- .smuce_test_partitions(n)
  valid <- vapply(parts, function(ends) {
    starts <- c(1L, head(ends, -1L) + 1L)
    all(vapply(seq_along(ends), function(i) {
      is.finite(smuce_cost(y[starts[i]:ends[i]], gamma, sigma2, n))
    }, logical(1)))
  }, logical(1))
  segment_counts <- vapply(parts, length, integer(1))
  min_segments <- min(segment_counts[valid])
  candidates <- parts[valid & segment_counts == min_segments]
  costs <- vapply(candidates, .smuce_test_partition_cost, numeric(1),
    y = y, gamma = gamma, sigma2 = sigma2, n = n
  )

  actual <- svp_smuce(y, gamma, sigma2)
  actual_cost <- .smuce_test_partition_cost(
    y, actual, gamma, sigma2, n
  )

  expect_equal(length(actual), min_segments)
  expect_equal(actual_cost, min(costs), tolerance = 1e-10)
  expect_equal(tail(actual, 1), n)
})

test_that("R and C++ SMUCE agree on generated series of length 100", {
  for (sigma2 in c(0.5, 1, 2)) {
    set.seed(203)
    y <- ts_generator(
      chpts = c(30L, 65L, 100L),
      parameters = c(0, 2.5, -1),
      sd_noise = sqrt(sigma2),
      type = "gauss"
    )
    expect_length(y, 100L)

    for (gamma in c(0, 0.75, 1.5, 3)) {
      expect_identical(
        svp_smuce(y, gamma, sigma2),
        svp_smuce_cpp(y, gamma, sigma2)
      )
    }
  }
})

test_that("SMUCE is invariant under additive shifts", {
  set.seed(20260912)
  y <- c(rnorm(30), 2 + rnorm(30), -1 + rnorm(30))
  gamma <- 1
  sigma2 <- 1

  expected_r <- svp_smuce(y, gamma, sigma2)
  expected_cpp <- svp_smuce_cpp(y, gamma, sigma2)
  expect_identical(expected_r, expected_cpp)

  # These shifts are large enough to exercise prefix-sum cancellation while
  # remaining exactly representable in the translated observations.
  for (offset in c(1e2, 1e4, 1e6)) {
    expect_identical(svp_smuce(y + offset, gamma, sigma2), expected_r)
    expect_identical(svp_smuce_cpp(y + offset, gamma, sigma2), expected_cpp)
    expect_equal(
      smuce_cost(y + offset, gamma, sigma2),
      smuce_cost(y, gamma, sigma2),
      tolerance = 1e-10
    )
  }
})

test_that("C++ SMUCE agrees with stepFit on generated series of length 1000", {
  skip_if_not_installed("stepR")

  cases <- list(
    constant = list(
      seed = 301L,
      chpts = 1000L,
      parameters = 0,
      sigma2 = 1,
      gamma = c(0, 1, 3)
    ),
    two_blocks = list(
      seed = 410L,
      chpts = c(500L, 1000L),
      parameters = c(0, 1),
      sigma2 = 1,
      gamma = c(0, 1.5, 3)
    ),
    three_blocks = list(
      seed = 411L,
      chpts = c(300L, 650L, 1000L),
      parameters = c(0, 1, -0.5),
      sigma2 = 1,
      gamma = c(0, 0.5, 1.5, 3)
    )
  )

  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    set.seed(case$seed)
    y <- ts_generator(
      chpts = case$chpts,
      parameters = case$parameters,
      sd_noise = sqrt(case$sigma2),
      type = "gauss"
    )
    expect_length(y, 1000L)

    for (gamma in case$gamma) {
      stepR_fit <- stepR::stepFit(
        y,
        q = gamma,
        family = "gauss",
        sd = sqrt(case$sigma2),
        intervalSystem = "all",
        lengths = seq_along(y),
        penalty = "sqrt"
      )
      expected <- as.integer(stepR_fit$rightIndex)

      expect_identical(
        svp_smuce_cpp(y, gamma, case$sigma2),
        expected,
        info = paste(case_name, "gamma =", gamma)
      )
    }
  }
})
