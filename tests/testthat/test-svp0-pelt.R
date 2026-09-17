run_pelt_pair <- function(data, gamma, test, subtests = "both") {
  list(
    no_pelt = svp0(
      data,
      gamma = gamma,
      test = test,
      subtests = subtests,
      PELT_pruning = FALSE
    ),
    pelt = svp0(
      data,
      gamma = gamma,
      test = test,
      subtests = subtests,
      PELT_pruning = TRUE
    )
  )
}

test_that("PELT preserves a block partition in both mode", {
  data <- c(rep(0, 5), rep(3, 5))
  range_valid <- function(segment, gamma) {
    diff(range(segment)) <= gamma
  }
  result <- run_pelt_pair(data, gamma = 0, test = range_valid)
  no_pelt <- result$no_pelt
  pelt <- result$pelt

  expect_equal(no_pelt$changepoints, c(5, 10))
  expect_equal(pelt$changepoints, no_pelt$changepoints)
  expect_equal(as.numeric(no_pelt$nb), c(1, 2, 3, 4, 5, 6, 2, 3, 4, 5))
  expect_equal(as.numeric(pelt$nb), as.numeric(no_pelt$nb))
  expect_equal(length(pelt$lastIndexSet), length(no_pelt$lastIndexSet))
})

test_that("both-mode subtests reduce candidates, not PELT", {
  data <- c(
    0.269605982037527392, -0.629985414039274993,
    0.868659827652264904, 1.727195517115240886,
    0.024187641767722562, 0.368025176992510361,
    -1.309204298208527018, 0.738621930711406027,
    0.044872987374990611, -1.048397200121202077,
    1.727851090445738036, -1.178599739022478721,
    0.653206711096908643, -0.368566491703636978,
    -0.599554643544779053
  )
  quantile_valid <- function(segment, gamma) {
    valid_QUANTILE(segment, gamma)
  }
  result <- run_pelt_pair(data, gamma = 2, test = quantile_valid)
  no_pelt <- result$no_pelt
  pelt <- result$pelt
  no_subtests <- svp0(
    data,
    gamma = 2,
    test = quantile_valid,
    subtests = "none",
    PELT_pruning = FALSE
  )

  expect_equal(no_pelt$changepoints, c(2, 4, 10, 11, 15))
  expect_equal(pelt$changepoints, no_pelt$changepoints)
  expect_equal(
    as.numeric(no_pelt$nb),
    c(1, 2, 3, 4, 3, 4, 5, 4, 5, 6, 7, 2, 2, 3, 4)
  )
  expect_equal(as.numeric(pelt$nb), as.numeric(no_pelt$nb))
  expect_lt(sum(no_pelt$nb), sum(no_subtests$nb))
  expect_equal(length(no_pelt$lastIndexSet), 5L)
  expect_equal(length(pelt$lastIndexSet), length(no_pelt$lastIndexSet))
  expect_lt(length(no_pelt$lastIndexSet), length(no_subtests$lastIndexSet))
})

test_that("PELT preserves a selected generic partition in both mode", {
  data <- c(
    -1.12716761108375851, 0.96260949243829663,
    1.40646401623188866, -1.64164950889970007,
    -1.12690419254923113, 0.59154548050242317,
    -0.95630911950593578, -0.69152988634639334
  )
  generic_test <- function(segment, gamma) {
    length(segment) %% 3 == 0 || sum(segment) > 0
  }
  result <- run_pelt_pair(data, gamma = 0, test = generic_test)
  no_pelt <- result$no_pelt
  pelt <- result$pelt

  expect_equal(no_pelt$changepoints, c(1, 4, 5, 6, 7, 8))
  expect_equal(pelt$changepoints, no_pelt$changepoints)
  expect_equal(as.numeric(no_pelt$nb), c(1, 2, 2, 3, 2, 2, 2, 2))
  expect_equal(as.numeric(pelt$nb), as.numeric(no_pelt$nb))
  expect_equal(length(pelt$lastIndexSet), length(no_pelt$lastIndexSet))
})
