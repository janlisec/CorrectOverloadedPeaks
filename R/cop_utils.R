#' @title verify_suggested.
#' @description Check if packages are available and stop function otherwise.
#' @param pkg Package names to be checked.
#' @return NULL.
#' @keywords internal
#' @noRd
verify_suggested <- function(pkg) {
  # verify that suggested packages are available
  check_pkg <- sapply(pkg, requireNamespace, quietly = TRUE)
  if (!all(check_pkg)) {
    msg <- paste0(
      "The use of this function requires package", ifelse(sum(!check_pkg)>1, "s", ""),
      paste(names(check_pkg)[!check_pkg], collapse=", "),
      ". Please install."
    )
    stop(msg)
  }
  invisible(NULL)
}

#' @title is.FlatTopPeak.
#' @description Check if packages are available and stop function otherwise.
#' @param x Numeric vector, ion intensities of a peak to test for going into saturation.
#' @param thr Defines the fraction of data points in a peak which are higher than 0.9*max(peak) --> for a nice Gaussian peak definitely <0.2.
#' @return TRUE or FALSE.
#' @keywords internal
#' @noRd
is.FlatTopPeak <- function(x, thr=0.2) {
  k <- table(x > 0.9*max(x))
  return((sum(k)/length(k))>thr)
}

#' @title GroupByGaps.
#' @description .
#' @param x Numeric vector.
#' @param gap Difference up from which a new time group is assumed.
#' @return A factor vector of length(x) providing levels of x falling in the 
#'     same bin, where bins are determined by gap size.
#' @examples
#' x <- c(0, 1, 2, 0.1, 2.1, 0.5)
#' GroupByGaps(x, 0.3)
#' split(x, GroupByGaps(x, 0.3))
#' split(x, GroupByGaps(x, 0.4))
#' @keywords internal
#' @noRd
GroupByGaps <- function(x, gap) {
  stopifnot(is.numeric(x))
  idx <- rank(x)
  x <- x[order(x)]
  x <- c(T, diff(x)>gap)
  x <- factor(rep(1:sum(x), times=diff(c(which(x),length(x)+1))))
  return(x[idx])
}