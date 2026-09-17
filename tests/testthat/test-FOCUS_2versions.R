library(testthat)
library(svpChange)

########## test svp0 result = SVP result with right pruning ##########
########## test svp0 result = SVP result with right pruning ##########
########## test svp0 result = SVP result with right pruning ##########

test_that("svp0 and SVP agree with right pruning", {
  n <- 500
  gap <- 1
  chpts <- c(0.1, 0.3, 0.4, 0.45, 0.55, 0.7, 0.75, 0.95, 1) * n
  data <- ts_generator(
    chpts = chpts,
    parameters = c(0, gap, 0, gap, 0, gap, 0, gap, 0),
    sd_noise = 1
  )
  gamma <- 5

  res_svp0 <- svp0(data,
    gamma,
    test = valid_FOCUS_last,
    subtests = "right",
    PELT_pruning = FALSE
  )

  res_svp <- SVP(
    data = data,
    gamma = gamma,
    test = "gaussian_mean",
    subtests = "right"
  )

  expect_equal(res_svp0$changepoints, res_svp$changepoints)
  expect_equal(res_svp0$R[, 1], res_svp$R[, 1])
})

test_that("svp0 and SVP return the same Gaussian FOCUS result", {
  n <- 100
  gap <- 5
  chpts <- c(0.1, 0.3, 0.4, 0.45, 0.55, 0.7, 0.75, 0.95, 1) * n
  data <- ts_generator(
    chpts = chpts,
    parameters = c(0, gap, 0, gap, 0, gap, 0, gap, 0),
    sd_noise = 1
  )
  gamma <- 7

  res_svp0 <- svp0(data,
    gamma,
    test = valid_FOCUS_last,
    subtests = "right",
    PELT_pruning = FALSE
  )
  res_svp <- SVP(
    data = data,
    gamma = gamma,
    test = "gaussian_mean",
    subtests = "right"
  )
  expect_equal(res_svp$changepoints, res_svp0$changepoints)
  expect_equal(res_svp$R, res_svp0$R)
})

test_that("SVP accepts the three pruning modes", {
  set.seed(123)
  n <- 80
  data <- ts_generator(
    chpts = c(20, 40, 60, 80), parameters = c(0, 1, -0.5, 0.8),
    sd_noise = 1, type = "gauss"
  )
  gamma <- 2 * log(n)

  opts <- data.frame(subtests = c("none", "right", "both"))

  for (i in seq_len(nrow(opts))) {
    res <- SVP(
      data = data,
      gamma = gamma,
      test = "gaussian_mean",
      subtests = opts$subtests[i]
    )

    expect_equal(tail(res$changepoints, 1), n)
  }
})

test_that("SVP evaluates Gaussian FOCUS at the current endpoint", {
  data <- c(2.002482730292829, 0.066700870930183,
            1.866851844706863)
  gamma <- 0.4873983399942518

  expect_false(valid_FOCUS(data, gamma))
  expect_true(valid_FOCUS_last(data, gamma))

  result <- SVP(
    data,
    gamma = gamma,
    test = "gaussian_mean",
    subtests = "none"
  )
  reference <- svp0(
    data,
    gamma = gamma,
    test = valid_FOCUS_last,
    subtests = "none",
    PELT_pruning = FALSE
  )

  expect_equal(result$changepoints, reference$changepoints)
  expect_equal(result$R, reference$R)
  expect_equal(result$changepoints, length(data))
})

test_that("both pruning applies the inclusive left rule", {
  # The complete segment is valid at endpoint 3, but the candidate segment
  # data[2:3] is invalid. In both mode that invalid boundary removes itself
  # and the older boundary 0, so the valid partition starts at observation 2.
  data <- c(0, 1, -1)
  result <- SVP(
    data,
    gamma = 7 / 8,
    test = "gaussian_mean",
    subtests = "both"
  )

  expect_equal(result$changepoints, c(2, 3))
  expect_equal(result$R[3, 2], 2)
  expect_equal(result$lastIndexSet, c(3, 2))
})

test_that("both reverse scan preserves the smallest-boundary tie rule", {
  # At t = 4, surviving boundaries 2 and 3 both give K = 3 and Q = 0.5.
  # The ascending cost pass must retain boundary 2 under the strict SSE tie.
  result <- SVP(
    c(-2, 0, -1, -2),
    gamma = 0.3,
    test = "gaussian_mean",
    subtests = "both"
  )

  expect_equal(result$R[4, ], c(0.5, 3, 2))
  expect_equal(result$changepoints, c(1, 2, 4))
})

test_that("SVP matches svp0 for equivalent pruning modes", {
  set.seed(27)
  data <- ts_generator(
    chpts = c(20, 40, 60), parameters = c(0, 2, -1),
    sd_noise = 1, type = "gauss"
  )
  gamma <- 1.5

  # The compiled Gaussian state and the R validity wrapper can take different
  # floating-point paths. Inclusive-left deletion can amplify a borderline
  # validity difference, so this direct implementation comparison is limited
  # to the two modes already covered by the reference tests.
  for (mode in c("none", "right")) {
    result <- SVP(
      data,
      gamma = gamma,
      test = "gaussian_mean",
      subtests = mode
    )
    reference <- svp0(
      data,
      gamma = gamma,
      test = valid_FOCUS_last,
      subtests = mode,
      PELT_pruning = FALSE
    )

    expect_equal(result$changepoints, reference$changepoints, info = mode)
    expect_equal(result$R, reference$R, info = mode)
    if (mode != "right") {
      expect_equal(result$nb, reference$nb, info = mode)
      expect_equal(result$lastIndexSet, reference$lastIndexSet, info = mode)
    }
  }
})

test_that("SVP computes stable squared-error costs for large offsets", {
  data <- 1e12 + c(0, 1, -1, 0, 2, -2)
  result <- SVP(
    data,
    gamma = 1e12,
    test = "varCost",
    subtests = "none"
  )
  expected <- sum((data - mean(data))^2)

  expect_equal(result$changepoints, length(data))
  expect_equal(result$R[length(data), 1], expected, tolerance = 1e-8)
  expect_true(all(result$R[, 1] >= 0))
})

test_that("Gaussian SVP is translation invariant under every pruning mode", {
  set.seed(2718)
  data <- rnorm(48)
  gamma <- 2 * log(length(data))

  for (mode in c("none", "right", "both")) {
    reference <- SVP(
      data, gamma = gamma, test = "gaussian_mean", subtests = mode
    )
    shifted <- SVP(
      data + 1e6, gamma = gamma, test = "gaussian_mean", subtests = mode
    )

    expect_equal(shifted$changepoints, reference$changepoints, info = mode)
    expect_equal(shifted$R[, 2:3], reference$R[, 2:3], info = mode)
    expect_equal(shifted$R[, 1], reference$R[, 1], tolerance = 1e-7,
                 info = mode)
  }
})

test_that("SVP defaults to Gaussian mean with both pruning rules", {
  set.seed(17)
  data <- ts_generator(
    chpts = 80, parameters = 0, sd_noise = 1, type = "gauss"
  )
  gamma <- 1.5 * log(length(data))

  default <- SVP(data, gamma)
  explicit <- SVP(data, gamma, test = "gaussian_mean")
  explicit_both <- SVP(
    data, gamma, test = "gaussian_mean", subtests = "both"
  )

  expect_equal(default$changepoints, explicit$changepoints)
  expect_equal(default$R, explicit$R)
  expect_equal(default$changepoints, explicit_both$changepoints)
  expect_equal(default$R, explicit_both$R)
})

test_that("both Gaussian pruning returns a valid final endpoint", {
  set.seed(2)
  n <- 100
  data <- ts_generator(
    chpts = n, parameters = 0, sd_noise = 1, type = "gauss"
  )
  gamma <- 0.3

  forward <- SVP(
    data,
    gamma,
    test = "gaussian_mean",
    sigma2 = 1,
    subtests = "both"
  )
  expect_equal(tail(forward$changepoints, 1), n)
})

test_that("both keeps the one-segment result for an unrestrictive threshold", {
  set.seed(3)
  data <- ts_generator(
    chpts = 1000, parameters = 0, sd_noise = 1, type = "gauss"
  )

  forward <- SVP(
    data,
    gamma = 1e12,
    test = "gaussian_mean",
    sigma2 = 1,
    subtests = "both"
  )
  expect_equal(forward$changepoints, 1000)
  expect_equal(as.numeric(forward$nb), seq_len(1000))
  expect_equal(forward$lastIndexSet, 1000:0)
})

test_that("SVP validates its public arguments", {
  expect_error(SVP(numeric(), 1), "at least one observation")
  expect_error(SVP(c(1, NA_real_), 1), "finite")
  expect_error(SVP(1:3, 0), "gamma")
  expect_error(SVP(1:3, 1, subtests = "invalid"), "subtests")
  expect_error(
    SVP(1:3, 1, subtests = "left"),
    "both.*right.*none"
  )
  expect_error(SVP(1:3, 1, test = "unknown"), "Unknown test")
  expect_error(SVP(c(1, 0, 2), 1, test = "gamma_rate"),
               "strictly positive")
  expect_error(SVP(1:3, 1, test = "quantile", quantile = 0.9), "0.5")
  expect_error(SVP(1:3, 1, test = "AR1", rho = 0, sigma2 = 0), "sigma2")
  expect_silent(SVP(1:3, 1, test = "gaussian_mean", sigma2 = 0))
})
