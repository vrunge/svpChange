\name{AR1_rho}
\alias{AR1_rho}
\title{Robust AR(1) Autocorrelation Estimate}
\usage{
AR1_rho(data)
}
\arguments{
  \item{data}{A numeric vector with at least three observations.}
}
\value{A stationary AR(1) coefficient in the interval \code{[-0.999, 0.999]}.}
\description{
Computes the robust estimator from equation (3) of Chakar et al. (2017), using
the squared ratio of the medians of absolute lag-two and lag-one differences,
minus one.
}
\references{
Chakar, S., Lebarbier, E., Levy-Leduc, C. and Robin, S. (2017).
A robust approach for estimating change-points in the mean of an AR(1)
process. \emph{Bernoulli}, 23(2), 1408--1447.
}
\examples{
set.seed(1)
data <- ts_generator(
  chpts = 60, parameters = 0, sd_noise = 1,
  rho = 0.6, type = "gaussAR1"
)
AR1_rho(data)
}
