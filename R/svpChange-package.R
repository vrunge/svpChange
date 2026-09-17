#' svpChange: Smallest Valid Partitioning for Change-Point Detection
#'
#' @description
#' Provides methods for detecting multiple change-points and segmenting
#' univariate time series using Smallest Valid Partitioning (SVP).
#'
#' SVP searches for a partition with the smallest number of contiguous
#' segments whose observations satisfy a user-defined or built-in validity
#' test. Among partitions with the same number of segments, it minimizes a
#' within-segment sum-of-squared-errors criterion.
#'
#' @details
#' The package provides two main SVP interfaces:
#'
#' \describe{
#'   \item{\code{\link{SVP}}}{A compiled implementation with built-in
#'   validity tests, including Gaussian mean, Gamma rate, variance,
#'   quantile, rank-based, and AR(1) tests.}
#'   \item{\code{\link{svp0}}}{A flexible implementation accepting an
#'   arbitrary R function of the form
#'   \code{function(segment, gamma)} as its validity test.}
#' }
#'
#' Validity-test helpers include:
#' \code{\link{valid_FOCUS}},
#' \code{\link{valid_AR1}},
#' \code{\link{valid_SSE}},
#' \code{\link{valid_RANGE}},
#' \code{\link{valid_RANGE_SLACK}},
#' \code{\link{valid_QUANTILE}},
#' \code{\link{valid_SCALE}}, and
#' \code{\link{valid_OP}}.
#'
#' The package also provides specialized SMUCE routines for fixed-variance
#' Gaussian observations. \code{\link{svp_smuce}} is an R implementation,
#' while \code{\link{svp_smuce_cpp}} provides a compiled implementation.
#' Both minimize the number of SMUCE-valid segments and then minimize the
#' constrained Gaussian residual cost.
#'
#' \code{\link{OP}}, \code{\link{PELT}}, and \code{\link{SN}} are available
#' for comparisons with other dynamic-programming algorithms.
#' \code{\link{AR1_rho}} and \code{\link{AR1_single_change}} provide AR(1)
#' diagnostics, and \code{\link{ts_generator}} generates simulation data.
#'
#' The computational implementations are written in C++ and exposed to R
#' through Rcpp. The R implementations of the SMUCE routines are retained
#' as transparent reference implementations.
#'
#' @references
#' Romano, G., Eckley, I. A., Fearnhead, P., and Rigaill, G. (2023).
#' Fast Online Changepoint Detection via Functional Pruning CUSUM Statistics.
#' \emph{Journal of Machine Learning Research}, 24(81), 1--36.
#' \url{https://www.jmlr.org/papers/v24/21-1230.html}
#'
#' Frick, K., Munk, A., and Sieling, H. (2014).
#' Multiscale Change-Point Inference.
#' \emph{Journal of the Royal Statistical Society: Series B}, 76(3),
#' 495--580.
#' \url{https://doi.org/10.1111/rssb.12047}
#'
#' Chakar, S., Lebarbier, E., Levy-Leduc, C., and Robin, S. (2017).
#' A robust approach for estimating change-points in the mean of an AR(1)
#' process.
#' \emph{Bernoulli}, 23(2), 1408--1447.
#'
#' @seealso
#' \code{\link{SVP}}, \code{\link{svp0}},
#' \code{\link{valid_FOCUS}},
#' \code{\link{svp_smuce}}
#'
#' @docType package
#' @name svpChange
#' @keywords package
"_PACKAGE"
