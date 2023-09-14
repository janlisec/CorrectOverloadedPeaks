is.FlatTopPeak <-
function(x, thr=0.2) {
  # x : numeric vector, ion intensities of a peak to test for going into saturation
  # thr : defines the fraction of datapoints in a peak which are higher than 0.9*max(peak) --> for a nice Gaussian peak definitly <0.2
  k <- table(x > 0.9*max(x))
  return((sum(k)/length(k))>thr)
}
