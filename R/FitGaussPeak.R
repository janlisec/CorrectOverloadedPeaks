#'@title Extrapolate a flat top peak using a Gauss approach.
#'
#'@description
#'\code{FitGaussPeak} will take retention time ('x') and intensity ('y') data and extrapolate all points above
#'a certain threshold based on further parameters using a Gaussian approach.
#'
#'@details
#'This function is mainly used internally (\code{\link{CorrectOverloadedPeaks}}) but can be of value on it's
#'own to test brute force peak reconstruction given that appropriate base peak chromatograms are available.
#'
#'@param x A numeric vector, retention times.
#'@param y A numeric vector, ion intensities.
#'@param scale_range Specifies the expected range for the true peak to exceed the observed, where scale_range=c(1,100) would assume anything between not overloaded and 100-fold overloaded.
#'@param steps Specifies a step parameter used to create a sequence within 'scale_range' to test for good fits, higher=more precision, fewer=faster.
#'@param cutoff Overloaded peaks will be screwed from Gaussian shape already when approaching detector saturation (DS), cutoff=0.95 ensures that points just before DS will not be used to model fit.
#'@param idx If not NULL, 'idx' is expected to specify points to correct explicitly (as a numeric-vector within 1:length(x)).
#'@param weight_front A weighting parameter to punish deviations in peak front and tail differently; 0.5=use front/tail equally, 1=use only front, 0=use only tail.
#'@param strip_data Use all provided data if 'none' (default). Strip 'front' or 'tail' data in case you observe peak fronting or tailing respectively.
#'@param account_for_baseline_offset If TRUE will subtract min(y) from y before fitting parameters.
#'@param method The method for peak shape calculation. Can be 'Gauss' or 'EMG' (exponentially modified gauss).
#'@param silent For testing purposes some QC-plot will be generated if silent=FALSE.
#'@param fix_sd Supply a fix standard deviation (sd) for the peak or leave NULL to estimate sd within function.
#'@param ... passed to the QC plot function, e.g. 'main' or 'xlab'.
#'
#'@return
#'An annotated plot of the mass spectrum and detailed information within the console (if silent=FALSE) and
#'the optimal fitted data points (vector of length(y), returned invisible).
#'
#'@examples
#'#load test data
#'data("mzXML_data")
#'names(mzXML_data)
#'str(mzXML_data[["scan"]][[1]])
#'pk <- ModelGaussPeak(height=10^7, width=3, scan_rate=10, e=0, ds=8*10^6, base_line=10^2)
#'plot(pk, main="Gaussian peak of true intensity 10^7 but cutt off at 8*10^6")
#'idx <- pk[,"int"]>0.005 * max(pk[,"int"])
#'tmp <- FitGaussPeak(x=pk[idx,"rt"], y=pk[idx,"int"], silent=FALSE, xlab="RT", ylab="Intensity")
#'
#'@export
#'
#'@importFrom grDevices grey
#'@importFrom graphics legend
#'@importFrom graphics lines
#'@importFrom graphics plot
#'@importFrom graphics points
#'@importFrom stats dnorm
#'@importFrom stats pnorm
#'@importFrom stats median
#'@importFrom stats sd
#'@importFrom stats optim
#'

FitGaussPeak <-
function(x, y, scale_range=c(1,10), steps=10, cutoff=0.95, idx=NULL, weight_front=0.5, strip_data="none", account_for_baseline_offset=TRUE, method=c("Gauss","EMG")[1], silent=TRUE, fix_sd=NULL, ...) {
  # potential furhter parameter
  # check_peak_shape [ToDo] could be included here but is done in CorrectOverloadedPeaks corrently.
  #helper function 'demg' : density of exponentially modified gauss
  demg <- function(x, mean=0, sd=1, l=1) {
    # erfc has to be defined as not standard in R
    erfc <- function(x) { 2 * pnorm(x * sqrt(2), lower.tail = FALSE) }
    return((l/2)*exp(l/2*(2*mean+l*sd^2-2*x))*erfc((mean+l*sd^2-x)/sqrt(2*sd)))
  }
  # set up optimization function, basically establish a parameter defined (exponetially modified) gauss peak 
  # and calculate difference from observed data (sum of squared error) weighted by w and by number of
  # corrected data points smaller than observed
  fn_g <- function(p) {
    d = p[3]*dnorm(x, mean=p[1], sd=p[2])
    e <- sqrt(sum(w*(d[-idx]-y[-idx])^2))/length(x[-idx])*(1+10*sum(x[idx]>d[idx]))
    #e <- sqrt(sum(w*(d[-idx]-y[-idx])^2))/length(x[-idx])*(1+10^sum(x[idx]>d[idx]))
    #e <- sqrt(sum(w*(c(0.7,0.3)[1+(d[-idx]>y[-idx])])*(d[-idx]-y[-idx])^2))/length(x[-idx])
    return(e)
  }
  fn_emg <- function(p) {
    d = p[3]*demg(x, mean=p[1], sd=p[2], l=p[4])
    e <- sqrt(sum(w*(d[-idx]-y[-idx])^2))/length(x[-idx])*(1+10*sum(x[idx]>d[idx]))
    return(e)
  }
  if (account_for_baseline_offset) {
    baseline_offset <- min(y)
    y <- y-baseline_offset
  } else {
    baseline_offset <- 0
  }
  y_max <- max(y)
  sr <- median(diff(x))
  # keep original x and y values (as they may be modified before fitting)
  x_original <- x; y_original <- y
  # filter for points relative to maximum or depending on idx
  if (is.null(idx)) {
    idx <- range(which(y>=(cutoff*y_max))); idx <- seq(idx[1],idx[2])
  }
  # modify idx in case that fronting or tailing is specified; keep unchanged data
  if (strip_data=="front") {
    # limit front values to the last one before overloading starts
    front_scans <- 1:ifelse(min(idx)<=1, 1, min(idx)-1)
    if (length(front_scans)>1) {
      flt <- front_scans[-length(front_scans)]
      x <- x[-flt]; y <- y[-flt]; idx <- idx-length(flt)
    }
  }
  # tailing
  if (strip_data=="tail") {
    # limit tail values to the last one before overloading starts
    tail_scans <- ifelse(max(idx)<length(x), max(idx)+1, length(x)):length(x)
    if (length(tail_scans)>1) {
      flt <- tail_scans[-1]
      x <- x[-flt]; y <- y[-flt]
    }
  }
  # weighting parameter for front/tail values
  w <- c(rep(weight_front,min(idx)-1),rep(1-weight_front,length(x)-max(idx)))
  # expected peak minimum as a boundary in 'optim'
  exp_min <- scale_range[1]*y_max
  exp_max <- scale_range[2]*y_max
  mu <- round(mean(x[idx]),1)
  # guess sigma from data or specify fix if strange values occur
  if (is.null(fix_sd)) {
    sig <- ifelse(abs(sd(x[idx])-1.01)>1 || length(idx)==1, 0.5, round(sd(x[idx]),2))
    #sig <- diff(range(x_original))/6 
  } else {
    sig <- fix_sd
  }
  if (length(idx)>=1) {
    # solution is sensitive to algorithm and starting parameters therefore scale is tested by brute force for anything within scale_range
    if (method=="EMG") test <- lapply(seq(exp_min,exp_max,length.out=steps), function(s) { stats::optim(par=c(mu,sig,s,l=1), fn=fn_emg, method="L-BFGS-B", lower=c(mu-0.5,0.01,exp_min,0.1), upper=c(mu+0.5,2,exp_max,10)) })
    if (method=="Gauss") test <- lapply(seq(exp_min,exp_max,length.out=steps), function(s) { stats::optim(par=c(mu,sig,s), fn=fn_g, method="L-BFGS-B", lower=c(mu-0.5,0.01,exp_min), upper=c(mu+0.5,2,exp_max)) })
    
    # keep only converging solutions
    #print(test)
    test_c <- sapply(test,function(x){x$convergence})==0
    if (any(test_c)) {
      # discard solutions which are non converging and x-fold worse than best sollution
      test <- test[which(test_c)]
      test_v <- sapply(test,function(x){x$value})
      flt <- which(test_v<(1.5*min(test_v)))
      test <- test[flt]
      # which is the best result (minimum error)
      best_res <- which.min(test_v[flt])
      # get optimized parameters, with minimal sum of squared error
      p.fit <- test[[best_res]]$par
      # get optimal fitted data
      if (method=="EMG") y_fit <- round(p.fit[3]*demg(x_original, mean=p.fit[1], sd=p.fit[2], l=p.fit[4]))
      if (method=="Gauss")  y_fit <- round(p.fit[3]*dnorm(x_original, mean=p.fit[1], sd=p.fit[2]))
      
      #print(test)
      if (!silent) {
        message(paste0("Number of converging sollutions: ", sum(test_c), ", keeping ", length(test)))
        if (method=="EMG") test_plot <- sapply(test, function(z) {p <- z$par; return(baseline_offset + p[3]*demg(x_original, mean=p[1], sd=p[2], l=p[4]))})
        if (method=="Gauss") test_plot <- sapply(test, function(z) {p <- z$par; return(baseline_offset + p[3]*dnorm(x_original, mean=p[1], sd=p[2]))})
        test_cols <- test_v[flt]/max(test_v[flt])
        if(min(test_cols)>0.1) test_cols <- grDevices::grey(test_cols-0.1)
        # establish plot
        plot(x=range(x_original), y=range(c(range(y_original), range(test_plot))), type="n", ...)
        # correct index in case front was removed
        if (strip_data=="front") { idx <- idx + length(x_original)-length(x) }
        # plot fitted lines
        for (i in 1:ncol(test_plot)) lines(x=x_original, y=test_plot[,i], col=test_cols[i])
        # plot best fit in green
        lines(x=x_original, y=test_plot[,best_res], col=3, lty=2, lwd=2)
        # plot original data without overloaded/corrected region
        points(x_original[-idx], baseline_offset+y_original[-idx], pch=21, bg=1, cex=1)
        # plot corrected data values
        points(x_original[idx], baseline_offset+y_original[idx], pch=21, bg=grDevices::grey(0.8), cex=1)
        # indicate tail and front values in case they were stripped
        if (strip_data=="tail") points(x_original[min(c((max(idx)+2),length(x_original))):length(x_original)], baseline_offset+y_original[min(c((max(idx)+2),length(x_original))):length(x_original)], pch=21, bg="white", cex=1)
        if (strip_data=="front") points(x_original[1:max(c(min(idx)-2),1)], baseline_offset+y_original[1:max(c(min(idx)-2),1)], pch=21, bg="white", cex=1)
        graphics::legend(x="topleft", legend=c(paste("method =",method), 
                                     paste("cutoff =", round(max(y_original[-idx]))),
                                     paste("baseline_offset =", round(baseline_offset)),
                                     paste("exp_min =", formatC(exp_min,2,format="e")),
                                     paste("exp_max =", formatC(exp_max,2,format="e")),
                                     paste("steps =", steps), 
                                     paste("mean_start =", mu),
                                     paste("sd_start =", round(sig,2))), bty="n", text.col=4)
        graphics::legend(x="topright", legend=c(paste("max_int =",formatC(max(test_plot[,best_res]),2,format="e")), 
                                      paste("err_fit_mean =",round(min(test[[best_res]]$value))), 
                                      paste("s_fit_start =",formatC(seq(exp_min,exp_max,length.out=steps)[which.min(test_v)],2,format="e")),
                                      paste("lambda_fit =",round(p.fit[4],1)),
                                      paste("mean_fit =",round(p.fit[1],1)),
                                      paste("sd_fit =",round(p.fit[2],2))), bty="n", text.col=4)
      }
      if (account_for_baseline_offset) {
        y_fit <- y_fit+baseline_offset
      }
      return(y_fit)
    } else {
      warning("[FitGaussPeak] No converging sollution found by 'optim'.")
      return(y_original)
    }
  } else {
    warning("[FitGaussPeak] No datapoints below cutoff*max(y).")
    return(y_original)    
  }
}
