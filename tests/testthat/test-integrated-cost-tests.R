test_that("SVP dispatches every former cost test", {
  set.seed(123)
  y <- rep(c(0, 1, -0.5), each = 20) + rnorm(60)
  tests <- c("quantileExact", "quantile", "varCost",
             "WilcoxonCost", "MedianMoodCost")
  for (test_name in tests) {
    result <- SVP(y, gamma = 10, test = test_name,
                  subtests = "right", quantile = 0.1)
    expect_true(is.list(result), info = test_name)
    expect_true(all(c("changepoints", "nb", "R") %in% names(result)),
                info = test_name)
    expect_equal(tail(result$changepoints, 1), length(y), info = test_name)
  }
})

test_that("quantile is only accepted on its valid parameter range", {
  y <- rnorm(30)
  expect_error(SVP(y, 5, "quantile", quantile = 0))
  expect_error(SVP(y, 5, "quantileExact", quantile = 1))
  expect_error(SVP(y, 5, "quantile", quantile = NA_real_))
  expect_silent(SVP(y, 5, "quantile", quantile = 0.1))
})

test_that("quantile tests evaluate short non-singleton segments", {
  y <- c(0, 10, 20)

  exact <- SVP(
    y, gamma = 5, test = "quantileExact",
    subtests = "none", quantile = 0.25
  )
  approximate <- SVP(
    y, gamma = 5, test = "quantile",
    subtests = "none", quantile = 0.25
  )

  expect_gt(exact$R[3, 2], 1)
  expect_gt(approximate$R[2, 2], 1)
})

test_that("varCost evaluates short non-singleton segments", {
  result <- SVP(c(0, 10), gamma = 1, test = "varCost",
                subtests = "none")

  expect_equal(result$changepoints, c(1, 2))
})

test_that("rank cost tests retain exact tie handling", {
  wilcoxon_valid <- function(x, gamma) {
    if (length(x) < 2) return(TRUE)
    best <- 0
    for (u in seq_len(length(x) - 1)) {
      score <- outer(
        x[seq_len(u)], x[(u + 1):length(x)],
        function(a, b) ifelse(a < b, 0.5, ifelse(a > b, -0.5, 0))
      )
      best <- max(best, abs(sum(score)))
    }
    best < gamma
  }

  mood_valid <- function(x, gamma) {
    n <- length(x)
    if (n < 2) return(TRUE)
    med <- sort(x, partial = n %/% 2 + 1)[n %/% 2 + 1]
    below <- x < med
    above <- x > med
    total_below <- sum(below)
    total_above <- sum(above)
    best <- 0
    for (u in seq_len(n - 1)) {
      a11 <- sum(below[seq_len(u)])
      a12 <- sum(above[seq_len(u)])
      a21 <- total_below - a11
      a22 <- total_above - a12
      nA <- a11 + a12
      nB <- a21 + a22
      if (nA > 0 && nB > 0 && total_below > 0 && total_above > 0) {
        determinant <- a11 * a22 - a12 * a21
        best <- max(
          best,
          (nA + nB) * determinant^2 /
            (nA * nB * total_below * total_above)
        )
      }
    }
    best < gamma
  }

  set.seed(1)
  data <- round(c(rnorm(20), rnorm(20, 1)), 1)
  for (gamma in c(3, 7, 15)) {
    for (subtests in c("none", "right")) {
      expect_equal(
        SVP(data, gamma, "WilcoxonCost", subtests = subtests)$changepoints,
        svp0(data, gamma, wilcoxon_valid,
             subtests = subtests, PELT_pruning = FALSE)$changepoints
      )
      expect_equal(
        SVP(data, gamma, "MedianMoodCost", subtests = subtests)$changepoints,
        svp0(data, gamma, mood_valid,
             subtests = subtests, PELT_pruning = FALSE)$changepoints
      )
    }
  }
})

test_that("integrated AR1 cost tests retain rho metadata", {
  set.seed(321)
  y <- numeric(80)
  mu <- rep(c(0, 1), each = 40)
  y[1] <- mu[1] + rnorm(1) / sqrt(1 - 0.6^2)
  for (i in 2:length(y))
    y[i] <- mu[i] + 0.6 * (y[i - 1] - mu[i - 1]) + rnorm(1)

  for (test_name in c("AR1", "AR1Profile", "AR1Focus")) {
    result <- SVP(y, 10, test_name, rho = 0.6, sigma2 = 1)
    expect_equal(result$rho, 0.6, info = test_name)
    expect_equal(result$sigma2, 1, info = test_name)
    expect_equal(tail(result$changepoints, 1), length(y), info = test_name)
  }
})

test_that("profiled AR1 tests ignore sigma2", {
  set.seed(322)
  y <- numeric(40)
  y[1] <- rnorm(1)
  for (i in 2:length(y)) y[i] <- 0.4 * y[i - 1] + rnorm(1)

  profile_default <- SVP(
    y, gamma = 6, test = "AR1Profile", rho = 0.4, sigma2 = 1
  )
  profile_zero <- SVP(
    y, gamma = 6, test = "AR1Profile", rho = 0.4, sigma2 = 0
  )
  profile_flag <- SVP(
    y, gamma = 6, test = "AR1", rho = 0.4, sigma2 = 0,
    profile_sigma = TRUE
  )

  expect_equal(profile_zero$changepoints, profile_default$changepoints)
  expect_equal(profile_flag$changepoints, profile_default$changepoints)
  expect_error(
    SVP(y, gamma = 6, test = "AR1", rho = 0.4, sigma2 = 0),
    "sigma2"
  )
  expect_error(
    SVP(y, gamma = 6, test = "AR1Focus", rho = 0.4, sigma2 = 0),
    "sigma2"
  )
})
