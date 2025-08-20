

devtools::install_github("alexandre-combeau/changepoints")
livrary(changepoints)
library(svpChange)

n <- 1000
gap <- 10
chpts = c(0.3,0.5,0.7,1)*n
data <- tsGenerator(chpts = chpts,
                    parameters = c(0,gap,0,gap),
                    sdNoise = 0.01)
plot(data)
penalty <- 2 * log(n)

### OP
OPres <- OP(data, penalty)
OPres$changepoints

### SVP0 valid_RANGE
### SVP0 valid_RANGE
### SVP0 valid_RANGE
plot(data)
res_svp0 <- svp0(data, gamma = 4, test = valid_RANGE)
res_svp0$changepoints
res_svp0$nb

res_svp <- smallest_valid_partitioning_rcpp(data = data,
                                            gamma = 20,
                                            test = valid_RANGE)

res_svp$nb

res_svp$changepoints

data <- rep(c(0,1,0,1,2,0), each = 100)
plot(data)
res_svp0 <- svp0(data, gamma = 0.999, test = valid_RANGE)
res_svp0$changepoints
res_svp0 <- svp0(data, gamma = 1, test = valid_RANGE)
res_svp0$changepoints
res_svp0 <- svp0(data, gamma = 1.9999, test = valid_RANGE)
res_svp0$changepoints

res_svp0 <- svp0(data, gamma = 0.5, test = valid_RANGE, prune_if_PELT = TRUE)
res_svp0$changepoints
res_svp0$nb



res_svp0 <- svp0(data,
                 gamma = 5,
                 test = valid_FOCUS,
                 prune_if_unvalid = F,
                 prune_if_PELT = T)
res_svp0$changepoints
res_svp0$nb


res_svp <- smallest_valid_partitioning_rcpp(data = data,
                                            gamma = 20,
                                            test = valid_RANGE)

res_svp$changepoints











  n <- 500
  gap <- 2
  chpts = c(0.1,0.3,0.4,0.45,0.55,0.7,0.75,0.95,1)*n
  data <- tsGenerator(chpts = chpts,
                           parameters = c(0,gap,0,gap,0,gap,0,gap,0),
                           sdNoise = 1.5)
  gamma <- 10

  # Use valid_FOCUS as the test function
  svp_1 <- smallest_valid_partitioning(data, gamma, test = valid_FOCUS)
  svp_1$changepoints
  svp_2 <- smallest_valid_partitioning_rcpp(data, gamma, test = valid_FOCUS)
  svp_2$changepoints

  svp_3 <- svp0(data, gamma, test = valid_FOCUS)
  svp_3$changepoints






################################################################################

n <- 500
gap <- 3
chpts = c(0.1,0.3,0.4,0.45,0.55,0.7,0.75,0.95,1)*n
data <- tsGenerator(chpts = chpts,
                      parameters = c(0,gap,0,gap,0,gap,0,gap,0),
                      sdNoise = 1)
gamma <- 3
plot(data)

svp <- svp0(data, gamma, test = valid_QUANTILE,
            prune_if_unvalid = F,
            prune_if_PELT = T)
abline(v = svp$changepoints)
svp$nb

plot(svp$R[,2], type = 'l')
plot(svp$nb, type = 'l')
svp$R[svp$lastIndexSet,2]



