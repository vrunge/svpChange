
library(svpChange)
library(dust)

n <- 1000
gap <- 1
chpts = c(0.1,0.3,0.4,0.45,0.55,0.7,0.75,0.95,1)*n
data <- dataGenerator_1D(chpts = chpts,
                         parameters = c(0,gap,0,gap,0,gap,0,gap,0),
                         sdNoise = 1)
plot(data)
penalty <- 2 * log(length(data))

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


