focus_lil_reference <- function(y, sigma2 = 1) {
  n <- length(y)
  best <- -Inf
  for (a in 0:(n - 2L)) {
    for (c in (a + 2L):n) {
      m <- c - a
      local <- -Inf
      largest_penalty <- -Inf
      for (b in (a + 1L):(c - 1L)) {
        left <- b - a
        right <- c - b
        h <- max(min(left, m / left), min(right, m / right))
        z <- (mean(y[(a + 1L):b]) - mean(y[(b + 1L):c])) *
          sqrt(left * right / (m * sigma2))
        local <- max(local, z^2)
        largest_penalty <- max(largest_penalty,
                               8 * log(log(exp(1) * h)))
      }
      best <- max(best, local - largest_penalty)
    }
  }
  best
}

test_that("FOCUS LIL agrees with the published split penalty", {
  set.seed(184)
  for (i in seq_len(15L)) {
    y <- rnorm(8L)
    critical <- focus_lil_reference(y)
    expect_length(SVP(y, critical + 1e-7,
                      subtests = "both_focus_lil")$changepoints, 1L)
    expect_gt(length(SVP(y, critical - 1e-7,
                         subtests = "both_focus_lil")$changepoints), 1L)
  }
})

test_that("FOCUS LIL handles variance, location, and input validation", {
  y <- c(-1.1, 0.3, 1.5, 0.4, -0.8, 0.9, 1.3, -0.2)
  q <- focus_lil_reference(y) - 0.1
  expected <- SVP(y, q, subtests = "both_focus_lil")$changepoints
  expect_equal(SVP(y + 1e5, q,
                   subtests = "both_focus_lil")$changepoints, expected)
  expect_equal(SVP(2 * y, q, subtests = "both_focus_lil",
                   sigma2 = 4)$changepoints, expected)
  expect_equal(SVP(3, -10, subtests = "both_focus_lil")$changepoints, 1)
  expect_error(SVP(y, Inf, subtests = "both_focus_lil"), "finite")
  expect_error(SVP(y, q, subtests = "both_focus_lil", sigma2 = 0),
               "positive")
  expect_error(SVP(y, q, test = "gamma_rate",
                   subtests = "both_focus_lil"), "gaussian_mean")
})
