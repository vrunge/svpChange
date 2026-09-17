\name{valid_AR1}
\alias{valid_AR1}
\title{AR(1) Mean-Change Validity Test}
\usage{
valid_AR1(
  y,
  gamma,
  rho = NA_real_,
  sigma2 = 1,
  profile_sigma = FALSE
)
}
\arguments{
  \item{y}{A numeric AR(1) segment. At least four observations are needed for
  a non-trivial single-change scan.}
  \item{gamma}{A threshold for the AR(1) likelihood-ratio statistic.}
  \item{rho}{An optional known AR(1) coefficient. If \code{NA}, it is
  estimated robustly from \code{y}.}
  \item{sigma2}{A positive innovation variance used by the known-variance
  statistic. In the AR(1) model, this is the variance of the new random shock
  after accounting for the previous observation and the AR(1) mean structure;
  it is not the marginal variance of the observed series. It is ignored when
  \code{profile_sigma = TRUE}.}
  \item{profile_sigma}{If \code{TRUE}, estimate the innovation variance
  separately under the no-change and one-change models from their residual sums
  of squares. This is useful when the innovation scale is unknown.}
}
\value{
\code{TRUE} if the AR(1) statistic is strictly below \code{gamma}; otherwise,
\code{FALSE}.
}
\description{
Tests whether a segment is valid under a conditional Gaussian AR(1)
single-mean-change statistic implemented in R.
}
\details{
The function conditions on the first observation, transforms the remaining
observations into AR(1) innovations, and scans change locations that leave at
least two observations on each side. The
statistic is one half of the Gaussian likelihood-ratio statistic, using the
same convention as \code{\link{AR1_single_change}}. If \code{rho} is
\code{NA}, the robust estimator requires at least three observations. The
function can be passed to \code{\link{svp0}} directly when \code{rho} is fixed,
or through a closure when other optional arguments are needed. The caller
should enforce any minimum segment length required by the chosen SVP procedure.
When \code{profile_sigma = FALSE}, the statistic uses the supplied
\code{sigma2}; when \code{profile_sigma = TRUE}, the innovation variance is
profiled out and estimated under each competing model.
}
\examples{
set.seed(1)
data <- ts_generator(
  chpts = 60, parameters = 0, sd_noise = 1,
  rho = 0.6, type = "gaussAR1"
)
valid_AR1(data, gamma = 10, rho = 0.6, sigma2 = 1)
}
\seealso{
\code{\link{AR1_single_change}}, \code{\link{svp0}}
}
