
library(svpChange)
library(dust)

n <- 1000
gap <- 1
chpts = c(0.1,0.3,0.4,0.45,0.55,0.7,0.75,0.95,1)*n
data <- tsGenerator(chpts = chpts,
                    parameters = c(0,gap,0,gap,0,gap,0,gap,0),
                    sdNoise = 1)
plot(data)
penalty <- 2 * log(n)

### OP
OPres <- OP(data, penalty)
OPres$changepoints

### PELT
PELTres <- PELT(data, penalty)
PELTres$changepoints

### DUST
DUSTres <- dust.1D(data, penalty = penalty/2)
DUSTres$changepoints


### La bonne comparaison des coûts pour obtenir 0
2*(DUSTres$costQ + cumsum(data^2)/2) - OPres$costQ
PELTres$costQ  - OPres$costQ

plot(2*(DUSTres$costQ + cumsum(data^2)/2) - OPres$costQ)




################################################################################
####
#### SIMULATIONS
####


library(changepoints)
library(microbenchmark)


n <-  10^5
penalty <- 2 * log(n)/3
data <- tsGenerator(n, sdNoise = 1)
res <- PELT(data, penalty)
res$changepoints

####
#### MORE SIMULATIONS
####


library(changepoints)
library(svpChange)
library(microbenchmark)

n <-  10^5
penalty <- 2 * log(n)/3
nrep <- 10L

# Pre-generate data
datasets <- replicate(nrep, tsGenerator(n, sdNoise = 1), simplify = FALSE)

# External counter
counter1 <- 0
counter2 <- 0

bench <- microbenchmark(
  PELT_Runge = {
    counter1 <<- counter1 + 1
    PELT(datasets[[counter1]], penalty)
  },
  PELT_Rcpp = {
    counter2 <<- counter2 + 1
    pelt_rcpp(datasets[[counter2]], penalty)
  },
  times = nrep
)

print(bench)
