

library(testthat)
library(svpChange)

########## test PELT chgpt = OP chgpt ##########
########## test PELT chgpt = OP chgpt ##########
########## test PELT chgpt = OP chgpt ##########

test_that("test PELT result = OP result",
          {
            n <- 1000
            gap <- 1
            chpts = c(0.1,0.3,0.4,0.45,0.55,0.7,0.75,0.95,1)*n
            data <- tsGenerator(chpts = chpts,
                                parameters = c(0,gap,0,gap,0,gap,0,gap,0),
                                sdNoise = 1)
            penalty <- 2 * log(n)
            ### OP
            OPres <- OP(data, penalty)
            ### PELT
            PELTres <- PELT(data, penalty)

            expect_equal(OPres$changepoints, PELTres$changepoints)
          })



