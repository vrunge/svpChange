

#devtools::install_github("vrunge/svpChange")
library(svpChange)

# Function to test time complexity

test_time_complexity <- function(nb, gamma = 10)
{
  times <- numeric(length(nb))

  for (i in seq_along(nb))
  {
    print(nb[i])
    data <- as.vector(replicate(nb[i], c(rnorm(50), rnorm(50, 2))))
    times[i] <- system.time({res <- SVP(data, gamma, test = "gaussian_mean")})[[1]]
  }
  data.frame(nb = 100*nb, time = times)
}



# Define nb of pattern repetitions
nb <- 50:100

# Run the test
time_results <- test_time_complexity(nb)


###
### Plot results
###
library(ggplot2)
ggplot(time_results, aes(x = nb, y = time)) +
  geom_point() +
  geom_line() +
  scale_y_log10() +
  scale_x_log10() +
  labs(x = "Data length", y = "Elapsed time (s)",
       title = "Time complexity of SVP function")

###
### Compute log-log regression
###
log_nb <- log10(time_results$nb)
log_times <- log10(time_results$time)


###
### slope = 1 !!!
###
fit <- lm(log_times ~ log_nb)
slope <- coef(fit)[2]  # this is the slope
slope


