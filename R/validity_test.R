##  GPL-3 License
## Copyright (c) 2025 Alexandre Combeau

#' Validity Test Based on FOCuS
#'
#' @title Validity Test Based on FOCuS
#' @description Checks whether a given segment is valid (i.e., has no changepoint).
#' @param y A numeric vector representing a segment of the signal.
#' @param gamma A numeric value used as a threshold in the validation function and as a penalty for each segment.
#' @return A logical value indicating whether the segment is considered valid (TRUE if no changepoint detected).
#' @export
valid_FOCUS <- function(y, gamma)
{
  res <- FOCuS::FOCuS(y, gamma)

  if(res$changepoint == -1)
  {
    return(TRUE)
  }
  else
  {
    return(FALSE)
  }
}

################################################
#' Validity Test Based on FOCuS last statistic
#'
#' @title Validity Test Based on FOCuS last statistic
#' @description Checks whether a given segment is valid (i.e., has no changepoint).
#' @param y A numeric vector representing a segment of the signal.
#' @param gamma A numeric value used as a threshold in the validation function and as a penalty for each segment.
#' @return A logical value indicating whether the segment is considered valid (TRUE if no changepoint detected).
#' @export
valid_FOCUS_last <- function(y, gamma)
{
  res <- FOCuS::FOCuS(y, gamma)

  if(res$changepoint == -1)
  {
    return(TRUE)
  }
  else
  {
    return(FALSE)
  }
}

################################################

#' Validity Test Based on Variance
#'
#' @title Validity Test Based on Variance
#' @description Checks if the variance of a segment is lower than or equal to a threshold.
#' @param y A numeric vector representing a segment of the signal.
#' @param gamma A numeric threshold for the maximum allowed variance.
#' @return TRUE if the variance is less than or equal to gamma, or if the segment is very short.
#' @export
valid_VAR <- function(y, gamma)
{
  if(length(y) < 6)
  {
    return(TRUE)
  }

  sum((y - mean(y))^2) <= gamma
}

################################################

#' Validity Test Based on Range
#'
#' @title Validity Test Based on Range
#' @description Checks whether the range (max - min) of the segment is below a threshold.
#' @param y A numeric vector representing a segment of the signal.
#' @param gamma A numeric threshold for the maximum allowed range.
#' @return TRUE if the range is less than or equal to gamma.
#' @export
valid_RANGE <- function(y, gamma)
{
  max(y) - min(y) <= gamma
}

################################################

#' Validity Test Based on Trimmed Range
#'
#' @title Validity Test Based on Trimmed Range
#' @description Applies a slack range test using trimmed minimum and maximum (ignores smallest and largest 3 values).
#' @param y A numeric vector representing a segment of the signal.
#' @param gamma A numeric threshold for the trimmed range.
#' @return TRUE if the trimmed range is less than or equal to gamma, or if the segment is too short.
#' @export
valid_RANGE_SLACK <- function(y, gamma)
{
  if(length(y) < 7)
  {
    return(TRUE)
  }
  sortedY <- sort(y)
  sortedY[length(y) - 3] - sortedY[3] <= gamma
}

################################################

#' Validity Test Based on Interquantile Range
#'
#' @title Validity Test Based on Interquantile Range
#' @description Checks whether the interquantile range (95% - 5%) is below a threshold.
#' @param y A numeric vector representing a segment of the signal.
#' @param gamma A numeric threshold for the interquantile range.
#' @return TRUE if the interquantile range is less than or equal to gamma.
#' @export
valid_QUANTILE <- function(y, gamma)
{
  unname(diff(quantile(y, probs = c(0.05, 0.95)))) <= gamma
}

################################################
#' Optimal Partitioning Cost Test
#'
#' @title Optimal Partitioning Cost Test
#' @description Tests whether the total cost of the segment is smaller than the best two-part penalized segmentation cost.
#' @param y A numeric vector representing a segment of the signal.
#' @param gamma A numeric threshold for penalizing the introduction of a new segment
#' @return TRUE if the segment cannot be split into two parts with lower (penalized) cost.
#' @export
valid_OP <- function(y, gamma)
{
  len <- length(y)
  if(len == 1){return(TRUE)}
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
  return(test = (total < (val + gamma)))
}



