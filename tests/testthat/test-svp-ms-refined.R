refined_critical_reference <- function(y, sigma2 = 1) {
  n <- length(y)
  best <- -Inf
  for (a in 0:(n - 2L)) {
    for (c in (a + 2L):n) {
      m <- c - a
      for (b in (a + 1L):(c - 1L)) {
        left <- b - a
        right <- c - b
        z <- abs(mean(y[(a + 1L):b]) - mean(y[(b + 1L):c])) *
          sqrt(left * right / (m * sigma2))
        penalty <- sqrt(2 * (1 + log(m / min(left, right)) +
                               0.5 * log(n / m)))
        best <- max(best, z - penalty)
      }
    }
  }
  best
}

test_that("refined multiscale threshold matches the exact interval scan", {
  set.seed(381)
  for (i in seq_len(12L)) {
    y <- rnorm(8L)
    critical <- refined_critical_reference(y)
    expect_length(SVP(y, critical + 1e-7,
                      subtests = "both_ms_refined")$changepoints, 1L)
    expect_gt(length(SVP(y, critical - 1e-7,
                         subtests = "both_ms_refined")$changepoints), 1L)
  }
})

test_that("refined multiscale statistic uses the supplied variance", {
  y <- c(-1.1, 0.3, 1.5, 0.4, -0.8, 0.9, 1.3, -0.2)
  q <- refined_critical_reference(y) - 0.1
  expected <- SVP(y, q, subtests = "both_ms_refined")$changepoints
  expect_equal(SVP(y + 1e5, q, subtests = "both_ms_refined")$changepoints,
               expected)
  expect_equal(SVP(2 * y, q, subtests = "both_ms_refined",
                   sigma2 = 4)$changepoints, expected)
  expect_equal(SVP(3, -10, subtests = "both_ms_refined")$changepoints, 1)
  expect_error(SVP(y, Inf, subtests = "both_ms_refined"), "finite")
  expect_error(SVP(y, q, subtests = "both_ms_refined", sigma2 = 0),
               "positive")
  expect_error(SVP(y, q, test = "gamma_rate",
                   subtests = "both_ms_refined"), "gaussian_mean")
})
