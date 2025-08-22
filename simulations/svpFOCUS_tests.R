

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
plot(data)

res_svp0 <- svp0(data,
                 gamma,
                 test = valid_FOCUS,
                 prune_if_unvalid = T,
                 prune_if_PELT = T)
abline(v = res_svp0$changepoints)
res_svp0$nb
plot(res_svp0$nb, type = 'l')
res_svp0$R[res_svp0$lastIndexSet,2]


##########


res_svp <- svp(data = data,
               gamma = gamma,
               test = "gaussian_mean",
               prune_if_unvalid = F,
               prune_if_PELT = T)
abline(v = res_svp$changepoints)
res_svp$nb
plot(res_svp$nb, type = 'l')
res_svp$R[res_svp$lastIndexSet,2]
