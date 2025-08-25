

devtools::install_github("alexandre-combeau/changepoints")
library(changepoints)
library(svpChange)


################################################################################
##
## valid_FOCUS test
##

n <- 500
gap <- 1
chpts = c(0.1,0.3,0.4,0.45,0.55,0.7,0.75,0.95,1)*n
data <- tsGenerator(chpts = chpts,
                      parameters = c(0,gap,0,gap,0,gap,0,gap,0),
                      sdNoise = 1)
gamma <- 3
bool <- F
res_svp0 <- svp0(data,
                 gamma,
                 test = valid_FOCUS,
                 prune_if_unvalid = bool,
                 prune_if_PELT = F)

res_svp <- SVP(data = data,
               gamma = gamma,
               test = "gaussian_mean",
               prune_if_unvalid = bool)

res_svp0$changepoints
res_svp$changepoints

