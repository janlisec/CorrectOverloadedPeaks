#'@title Create and modify parameters of an artificial chromatographic peak.
#'
#'@description
#'\code{ModelGaussPeak} will create a potentially overloaded Gaussian peak of requested width and height.
#'
#'@details
#'The main task of \code{\link{ModelGaussPeak}} is to create peak data in Gaussian shape for testing. 
#'Width is meant in the chromatographic sense, i.e. the time between peak front and tail hitting the baseline.
#'
#'@param height True peak height (=intensity counts).
#'@param width Peak width in time units (preferably seconds).
#'@param scan_rate Is determining the resolution of data points per time unit (preferably seconds).
#'@param e Error term giving the percent amount of deviation from the ideal Gaussian curve for individual data points.
#'@param ds Detector saturation. Intensity values will be cut off at this point if requested.
#'@param base_line Defines if peak is supposed to have a higher base level.
#'
#'@return
#'Dataframe with columns 'rt' and 'int'.
#'
#'@examples
#'ylim <- c(0,10^7)
#'par(mfrow=c(1,5))
#'pk <- ModelGaussPeak(height=10^7, width=4, scan_rate=10, e=0, ds=10^7, base_line=10^2)
#'plot(pk,ylim=ylim,main="standard")
#'pk <- ModelGaussPeak(height=10^7, width=4, scan_rate=10, e=0, ds=8*10^6, base_line=10^2)
#'plot(pk,ylim=ylim,main="flat top")
#'pk <- ModelGaussPeak(height=10^7, width=4, scan_rate=10, e=0, ds=8*10^6, base_line=10^5)
#'plot(pk,ylim=ylim,main="high baseline")
#'pk <- ModelGaussPeak(height=10^7, width=4, scan_rate=10, e=0.05, ds=8*10^6, base_line=10^5)
#'plot(pk,ylim=ylim,main="e=5%")
#'pk <- ModelGaussPeak(height=10^7, width=4, scan_rate=5, e=0.05, ds=8*10^6, base_line=10^5)
#'plot(pk,ylim=ylim,main="sr=5")
#'
#'@export
#'
#'@importFrom graphics plot
#'@importFrom stats dnorm
#'@importFrom stats runif
#'

ModelGaussPeak <- function(height=10^7, width=4, scan_rate=10, e=0, ds=10^7, base_line=10^2) {
# time series depending on scan_rate and desired width
  x <- seq(-0.5*width, 0.5*width, 1/scan_rate)
  xborder <- min(which(height*dnorm(seq(-7,0,0.1))>base_line))
  #if (!is.finite(xborder)) browser()
  xborder <- seq(-7,0,0.1)[xborder]
# idealized data points according to dnorm
  y <- dnorm(seq(xborder,-1*xborder,length.out=length(x)))
# scaled to desired peak height
  y <- height*y/max(y)
# and modified by defined error
  if (e>0) y <- sapply(y, function(z) {z+runif(1,-e,e)*z})
# before detector saturation is applied
  y[y>ds] <- ds
  return(data.frame("rt"=x, "int"=y))
}