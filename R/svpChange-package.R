

#' svpChange: Smallest Valid Partitioning for multiple change-point detection
#'
#' @description
#' Tools for multiple change-point detection and segmentation of **univariate** sequences using the
#' *smallest valid partitioning* (SVP) algorithm.
#'
#' SVP builds a segmentation by enforcing a **local notion of segment validity**: a candidate segment is
#' retained only if it passes a user-specified validity test (for instance, a single change-point test).
#' If a change is detected, the segment is deemed **invalid** and must be split, preventing invalid
#' segments from appearing in the final partition.
#'
#' The computational core is implemented in C++ and exposed to R via **Rcpp** for performance.
#'
#' @details
#' SVP is designed for situations where you want **explicit, test-based control** of within-segment
#' homogeneity rather than relying solely on a global penalized objective.
#'
#'
#' @section Functions:
#' - \code{\link{svp0}}:
#' - More functions here.
#'
#' @docType package
#' @name svpChange
#' @keywords internal
"_PACKAGE"
