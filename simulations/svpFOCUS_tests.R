

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
gamma <- 20
bool <- T
res_svp0 <- svp0(data,
                 gamma,
                 test = valid_FOCUS,
                 prune_if_unvalid = bool)

res_svp <- SVP(data = data,
               gamma = gamma,
               test = "gaussian_mean",
               prune_if_unvalid = bool)
res_svp0$changepoints
res_svp$changepoints

res_PELT <- PELT(data = data,
             penalty = gamma)
res_PELT$changepoints



########################################################


####

n <- 500
data <- rnorm(n)
gamma <- 20

res <- FOCuS::FOCuS(data, gamma)
res
valid_FOCUS(data, gamma)

y <- data
len <- length(y)

total <- sum((y - mean(y))^2)
val <- Inf
for(i in 2:len)
{
  segment1 <- y[1:(i-1)]
  segment2 <- y[i:len]
  mean_seg1 <- mean(segment1)
  mean_seg2 <- mean(segment2)
  temp <- (sum((segment1 - mean_seg1)^2) + sum((segment2 - mean_seg2)^2))
  if(temp < val){val <- temp}
}
(total - val)/2


