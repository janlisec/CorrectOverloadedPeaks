#' @title read.mzData.
#'
#' @description
#' \code{read.mzData} will import mzData as xcmsRaw-class objects.
#'
#' @details
#' The main task of read.mzData functions is to import mzData files to R. 
#'     Currently `xcmsRaw` is supported as an output format. I created this 
#'     function for legacy reasons as the mzData import is no longer supported 
#'     by `mzR` and consequently `xcms` since 09/2021.
#'     This is a quick and dirty implementation. It will work only for mslevel=1
#'     and a fixed set of base64 encoding parameters (size = 4, endian = "big").
#'     However, feel free to send me an e-mail if you are interested in using
#'     the function but cant get it working.
#'
#' @param filename A mzData file as exported by `xcms::write.mzdata()`.
#' @param fmt Output format. Currently only xcmsRaw is supported.
#' @param verbose Print messages to console.
#' 
#' @return
#' A generic R object of class xcmsRaw.
#' 
#' @examples
#' \dontrun{
#'   data(mzXML_data)
#'   write.mzXML(mzXML = mzXML_data, filename = "test.mzXML")
#'   x <- xcms::xcmsRaw("test.mzXML",  profstep=0)
#'   xcms::write.mzdata(x, file="test.mzData")
#'   x2 <- read.mzData(filename = "test.mzData")
#'   identical(str(x), str(x2))
#'   identical(x@env$intensity, x2@env$intensity)
#'   identical(x@env$mz, x2@env$mz)
#'   identical(x@scanindex, x2@scanindex)
#'   file.remove(c('test.mzData', 'test.mzXML'))
#' } 
#'
#' @export
#' 
read.mzData <- function(filename, fmt = c("xcmsRaw"), verbose = FALSE) {
  
  fmt <- match.arg(fmt)
  
  if (!is.character(filename)) stop("read.mzData: Parameter 'filename' has to be a string.")
  if (length(filename) > 1) filename <- paste(filename, collapse = "") # combine characters into a string
  
  if (!file.exists(filename)) stop("read.mzData: File ", filename, " does not exist.")
  
  if (tolower(tools::file_ext(filename))!="mzdata") warning("read.mzData: File '", filename, "' does not have extension 'mzData'.")
  
  # verify that suggested packages are available
  verify_suggested(c("xcms", "xml2", "methods", "bitops"))

  xt <- xml2::read_xml(x = filename)
  if (is.null(xt)) { stop("read.mzData: xml2::read_xml returns NULL upon import of file '", filename, "'.") }
  
  xt <- xml2::as_list(xt)
  xt <- xt[["mzData"]][["spectrumList"]]
  
  if (length(xt)==0) { stop("read.mzData: file '", filename, "' seems not to contain a valid 'spectrumList' slot.") }
  
  if (!verbose) message("Processing file ", filename)
  
  # remove empty scans if present as xcmsRaw format can not handle these
  if (fmt %in% c("xcmsRaw")) {
    flt <- sapply(xt, function(xn) { length(xn$intenArrayBinary$data) })==1
    if (!any(flt)) { stop("read.mzData: file '", filename, "' contains only empty scans (filed $intenArrayBinary$data.") }
    xt <- xt[flt]
  }
  
  test_scan <- xt[[1]]$spectrumDesc$spectrumSettings$spectrumInstrument
  # check if time in seconds is exported
  rt_par <- sapply(test_scan, function(y) { any(unlist(attributes(y))=="PSI:1000039") })
  fac <- 0
  if (sum(rt_par)==1) {
    idx <- which(rt_par)
    fac <- 1
  } else {
    # check if time in minutes is exported
    rt_par <- sapply(test_scan, function(y) { any(unlist(attributes(y))=="PSI:1000038") })
    if (sum(rt_par)==1) {
      idx <- which(rt_par)
      fac <- 60
    }
  }
  if (fac==0) { stop("read.mzData: file '", filename, "' contains not RT information (PSI:1000038 or PSI:1000039).") }
  rt <- unname(sapply(xt, function(xn) {
    x <- xn$spectrumDesc$spectrumSettings$spectrumInstrument
    as.numeric(attributes(x[[idx]])$value)
  }))*fac

  # get parameters of base64 encoding probably used  
  test_data <- xt[[1]]$intenArrayBinary$data
  size.int <- switch(
    attr(test_data, "precision"),
    "32"=4,
    "64"=8,
    4
  )
  endian.int <- attr(test_data, "endian")
  # extract intensity vectors
  int <- unname(lapply(xt, function(xn) {
    base64decode(z = xn$intenArrayBinary$data[[1]], what = "numeric", size = size.int, signed = TRUE, endian = endian.int)
  }))
  
  # get parameters of base64 encoding probably used  
  test_data <- xt[[1]]$mzArrayBinary$data
  size.mz <- switch(
    attr(test_data, "precision"),
    "32"=4,
    "64"=8,
    4
  )
  endian.mz <- attr(test_data, "endian")
  # extract mz vectors
  mz <- unname(lapply(xt, function(xn) {
    base64decode(z = xn$mzArrayBinary$data[[1]], what = "numeric", size = size.mz, signed = TRUE, endian = endian.mz)
  }))
  
  pol_par <- sapply(test_scan, function(y) { any(unlist(attributes(y))=="PSI:1000037") })
  if (sum(pol_par)==1) {
    idx <- which(pol_par)
    pol <- unname(sapply(xt, function(xn) {
      x <- xn$spectrumDesc$spectrumSettings$spectrumInstrument
      attributes(x[[idx]])$value
    }))
  } else {
    message("read.mzData: file '", filename, "' contains not polarity information (PSI:1000037).")
    pol <- rep("unknown", length(xt))
  }

  out <- xcms::deepCopy(methods::new("xcmsRaw"))
  assign(x = "intensity", value = unlist(int, use.names = FALSE), envir = out@env)
  assign(x = "mz", value = unlist(mz, use.names = FALSE), envir = out@env)
  out@tic <- unname(sapply(int, sum, USE.NAMES = FALSE))
  out@scantime <- unname(rt)
  out@scanindex <- as.integer(c(0, cumsum(unname(sapply(int, length, USE.NAMES = FALSE))))[1:length(int)])
  out@polarity <- factor(tolower(pol), levels = c("negative", "positive", "unknown"))
  out@acquisitionNum <- as.integer(unname(sapply(xt, function(xn) { attr(xn, "id") })))
  out@mzrange <- range(out@env$mz)
  out@scanrange <- range(out@acquisitionNum)
  out@filepath <- xcms::xcmsSource(filename)
  
  return(out)
}
