\name{valid_SCALE}
\alias{valid_SCALE}
\title{Validity Test Based on Robust Scale}
\usage{
valid_SCALE(y, gamma)
}
\arguments{
  \item{y}{A numeric vector representing a segment of the signal.}
  \item{gamma}{A numeric threshold for the robust scale of the segment.}
}
\value{
\code{TRUE} if the normalized median absolute deviation is less than or equal
to \code{gamma}; otherwise, \code{FALSE}.
}
\description{
Checks whether the Gaussian-normalized median absolute deviation of a segment
is below a threshold.
}
\details{
The scale is computed with \code{\link[stats:mad]{stats::mad()}} using the
segment median as centre and the constant \code{1.4826}. This constant makes
the MAD estimate the standard deviation for Gaussian data. Because the test is
based on the complete segment, its result is not necessarily preserved when
observations are appended; use \code{subtests = "none"} and
\code{PELT_pruning = FALSE} in \code{\link{svp0}} unless the required
pruning properties have been established for the intended application.
}
\examples{
set.seed(1)
data <- ts_generator(
  chpts = 30, parameters = 0, sd_noise = 1, type = "gauss"
)
valid_SCALE(data, gamma = 2)
}
