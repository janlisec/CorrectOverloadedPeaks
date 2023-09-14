GroupByGaps <-
function(x, gap) {
  # x : numeric vector
  # gap : difference up from which a new time group is assumed
  stopifnot(is.numeric(x))
  idx <- rank(x)
  x <- x[order(x)]
  x <- c(T, diff(x)>gap)
  x <- factor(rep(1:sum(x), times=diff(c(which(x),length(x)+1))))
  return(x[idx])
}
