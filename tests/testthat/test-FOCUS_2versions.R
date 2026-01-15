

library(testthat)
library(svpChange)
library(FOCuS)

########## test svp0 result = SVP result with prune_if_unvalid == TRUE ##########
########## test svp0 result = SVP result with prune_if_unvalid == TRUE ##########
########## test svp0 result = SVP result with prune_if_unvalid == TRUE ##########

test_that("test svp0 result = SVP result with prune_if_unvalid == TRUE",
          {
            n <- 500
            gap <- 1
            chpts = c(0.1,0.3,0.4,0.45,0.55,0.7,0.75,0.95,1)*n
            data <- tsGenerator(chpts = chpts,
                                parameters = c(0,gap,0,gap,0,gap,0,gap,0),
                                sdNoise = 1)
            gamma <- 5
            bool <- TRUE

            res_svp0 <- svp0(data,
                             gamma,
                             test = valid_FOCUS, #valid_FOCUS_last
                             prune_if_unvalid = bool)

            res_svp <- SVP(data = data,
                           gamma = gamma,
                           test = "gaussian_mean",
                           prune_if_unvalid = bool)

            expect_equal(res_svp0$changepoints, res_svp$changepoints)
            expect_equal(res_svp0$R[,1], res_svp$R[,1])
          })


test_that("svp0 svp and svp_old all equal",
          {
            n <- 50
            gap <- 5
            chpts = c(0.1,0.3,0.4,0.45,0.55,0.7,0.75,0.95,1)*n
            data <- tsGenerator(chpts = chpts,
                                parameters = c(0,gap,0,gap,0,gap,0,gap,0),
                                sdNoise = 1)
            bool <- TRUE

            gamma <- 7

            plot(data)

            res_svp0 <- svp0(data,
                             gamma,
                             test = valid_FOCUS, #valid_FOCUS_last
                             prune_if_unvalid = bool)
            abline(v = res_svp0$changepoints, col = "red")

            res_svp <- SVP(data = data,
                           gamma = gamma,
                           test = "gaussian_mean",
                           prune_if_unvalid = bool)
            abline(v = res_svp$changepoints, col = "blue", lty = 2, lwd = 2)


            res_svp_old <- SVP_old(data = data,
                                   gamma = gamma,
                                   test = "gaussian_mean",
                                   prune_if_unvalid = bool)
            abline(v = res_svp_old$changepoints, col = "green", lty = 3)


            expect_equal(res_svp$R, res_svp_old$R)

          })






