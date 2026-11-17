test_that("constant dyadic thresholds match SVP both", {
  set.seed(929)
  series <- list(rnorm(30), c(rnorm(10), rnorm(10, 2), rnorm(10)),
                 rep(3, 16), c(0, 10, 0, rep(0, 5)))

  for (y in series) {
    bands <- max(1L, floor(log2(length(y) / 2)) + 1L)
    for (gamma in c(0.5, 2, 20)) {
      expected <- SVP(y, gamma, subtests = "both")
      actual <- SVP_both_dyadic(y, rep(gamma, bands))
      expect_identical(actual$changepoints, expected$changepoints)
      expect_identical(actual$lastIndexSet, expected$lastIndexSet)
      expect_identical(actual$nb, expected$nb)
      expect_equal(actual$R, expected$R)
      expect_null(actual$costQ)
    }
  }
})

test_that("dyadic thresholds are selected by candidate segment length", {
  y <- c(0, 10, rep(0, 6))
  # For n = 8 the bands are (4, 8], (2, 4], and (1, 2].
  permissive <- SVP_both_dyadic(y, c(100, 100, 100))
  strict_at_length_two <- SVP_both_dyadic(y, c(100, 100, 0.1))

  expect_length(permissive$changepoints, 1L)
  expect_gt(length(strict_at_length_two$changepoints), 1L)
  expect_equal(tail(permissive$changepoints, 1L), length(y))
  expect_equal(tail(strict_at_length_two$changepoints, 1L), length(y))
})

test_that("dyadic SVP validates data and gamma dimensions", {
  expect_error(SVP_both_dyadic(numeric(), 1), "at least one")
  expect_error(SVP_both_dyadic(c(0, NA), c(1, 1)), "finite")
  expect_error(SVP_both_dyadic(1:8, c(1, 1)), "one value per dyadic")
  expect_error(SVP_both_dyadic(1:8, c(1, 1, 0)), "finite and positive")
  expect_error(SVP_both_dyadic(1:8, c(1, Inf, 1)), "finite and positive")
  expect_equal(SVP_both_dyadic(3, 1)$changepoints, 1)
})

test_that("dyadic SVP is invariant to additive shifts", {
  set.seed(930)
  y <- rnorm(40)
  gamma <- rep(10, 5)
  expected <- SVP_both_dyadic(y, gamma)$changepoints

  for (offset in c(1e2, 1e4, 1e6)) {
    expect_identical(
      SVP_both_dyadic(y + offset, gamma)$changepoints,
      expected
    )
  }
})
