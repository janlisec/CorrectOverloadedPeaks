#'@title Extrapolate a flat top peak using isotopic ratios.
#'
#'@description
#'\code{FitPeakByIsotopicRatio} will take a data frame containing peak data for retention time ('RT'), as well as
#'mass and intensity information of M0, M+1 and M+2 and extrapolate all points above a certain threshold for Int_M0
#'based on further parameters using an IsotopicRatio approach.
#'
#'@details
#'Isotopic ratios within ion traces of molecules can be considered stable. If this ratio is changed because one molecule,
#'let's say the M+0, is exceeding the detector range while another (say M+1) is still quantifiable, we therefore may
#'attempt to modify M+0 by multiplying the values of M+1 with a constant (the stable ratio). This constant is determined
#'ideally from the values within the peak front.
#'As this function is mainly used internally (\code{\link{CorrectOverloadedPeaks}}), it is not very flexible with respect
#'to the input format. Please prepare a dataframe according to the parameter specifications or process a file using
#'\code{\link{CorrectOverloadedPeaks}} with testing=TRUE, which will generate a list structure of such dataframes.
#'
#'@param cor_df A data frame containing information about the overloaded area; columns=(Scan, RT, mz0, int0, mz1, int1, mz2, int2, modified).
#'@param idx If not NULL, 'idx' is expected to specify points to correct explicitly (as a numeric-vector within 1:length(x)).
#'@param silent For testing purposes some QC-plot will be generated if silent=FALSE.
#'
#'@return
#'An annotated plot of the mass spectrum and detailed information within the console.
#'Main result will be returned invisible.
#'
#'@export
#'
#'@importFrom graphics legend
#'@importFrom graphics lines
#'@importFrom graphics par
#'@importFrom graphics plot
#'@importFrom graphics points
#'@importFrom stats median
#'

FitPeakByIsotopicRatio <-
function(cor_df=NULL, idx=NULL, silent=TRUE){
  # potential parameters
  cor_fac_front=TRUE
  ensure_minimum=TRUE
  n_iso <- length(grep("mz",colnames(cor_df)))-1
  # apply iso-ratio correction based on front or tail value (!) -- tail value will be systematically lower by ~3-5% of total
  old_vals <- cor_df
  if (is.null(idx)) {
    idx <- range(which(cor_df[,"int0"]>=0.95*max(cor_df[,"int0"]))); idx <- seq(idx[1],idx[2])
  }
  lim <- min(cor_df[idx,"int0"])
  front_scans <- 1:(min(idx)-1)
  tail_scans <- (max(idx)+1):nrow(cor_df)
  # assuming that the last provided isotope is the first non-overloaded isotope
  if (cor_fac_front) cor_fac <- median(cor_df[front_scans,paste0("int",n_iso)]/cor_df[front_scans,"int0"],na.rm=T) else cor_fac <- median(cor_df[tail_scans,paste0("int",n_iso)]/cor_df[tail_scans,"int0"],na.rm=T)
  if (!is.finite(cor_fac)) cor_fac <- 1
  cor_df[idx,"int0"] <- round(cor_df[idx,paste0("int",n_iso)]/cor_fac)
  if (ensure_minimum) {
    for (i in 0:n_iso) {
      flt <- cor_df[,paste0("int",i)] < old_vals[,paste0("int",i)]
      if (any(flt)) cor_df[flt,paste0("int",i)] <- old_vals[flt,paste0("int",i)]
    }
  }
  # some plot output
  if (!silent) {
    plot(x=range(cor_df[,"RT"]), y=c(0, max(cor_df[,"int0"])), type="n", main=paste("mz =", round(median(cor_df[,"mz0"]),4)), xlab="RT", ylab="Int")
    for (i in rev(0:n_iso)) {
      lines(x=cor_df[,"RT"], y=cor_df[,paste0("int",i)], col=i+1)
      points(x=cor_df[,"RT"], y=old_vals[,paste0("int",i)], pch=rep(c(21,22,24,25),2)[i+1], bg=i+1)
    }
    points(x=cor_df[idx,"RT"], y=old_vals[idx,"int0"], pch=21, bg=grDevices::grey(0.9))
    graphics::legend(x="topleft", bty="n", text.col=4, legend=c(
      paste("method =","Isoratio"),
      paste("iso_used =",paste0("M+",n_iso)),
      paste(ifelse(cor_fac_front, "front_ratio =", "tail_ratio ="), round(100*cor_fac,2), "%")))
    graphics::legend(x="topright", legend=paste("max_int =", formatC(max(cor_df[,"int0"]), digits=2, format="e")), bty="n", text.col=4)
  }
  return(cor_df)
}