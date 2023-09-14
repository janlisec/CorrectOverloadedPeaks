#' @title base64decode.
#' @description `base64decode` is a copy of a similar function from the caTools 
#'     package as this package is about to be archived (07/2018).
#' @details `base64decode` will convert base64 encoded strings into R values.#'
#' @param z The base64 encoded string.
#' @param what Define output type of z (e.g. 'numeric').
#' @param size Encoding size (provide if you know it).
#' @param signed Parameter passed through to `readBin`.
#' @param endian Parameter passed through to `readBin`. 
#' @return Decoded value of z.#' 
#' @export
base64decode <- function(z, what, size = NA, signed = TRUE, endian = .Platform$endian) {
  
  verify_suggested(c("bitops"))
  
  # this is a copy from the caTools function as this package is about to be archived (07/2018)
  if (!is.character(z)) {
    stop("base64decode: Input argument 'z' is suppose to be a string")
  }
  if (length(z) == 1) {
    z <- strsplit(z, NULL)[[1]]
  }
  if (length(z) %% 4 != 0) {
    warning("In base64decode: Length of base64 data (z) not a multiple of 4.")
  }
  alpha <- "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/="
  alpha <- strsplit(alpha, NULL)[[1]]
  y <- match(z, alpha, nomatch = -1) - 1
  if (any(y == -1)) {
    stop("base64decode: Input string is not in Base64 format")
  }
  if (any(y == 64)) {
    y <- y[y != 64]
  }
  neByte <- length(y)
  nBlock <- ceiling(neByte / 4)
  ndByte <- 3 * nBlock
  if (neByte < 4 * nBlock) {
    y[(neByte + 1):(4 * nBlock)] <- 0
  }
  dim(y) <- c(4, nBlock)
  x <- matrix(as.integer(0), 3, nBlock)
  x[1, ] <- bitops::bitOr(bitops::bitShiftL(y[1, ], 2), bitops::bitShiftR(y[2, ], 4))
  x[2, ] <- bitops::bitOr(bitops::bitShiftL(y[2, ], 4), bitops::bitShiftR(y[3, ], 2))
  x[3, ] <- bitops::bitOr(bitops::bitShiftL(y[3, ], 6), y[4, ])
  x <- bitops::bitAnd(x, 255)
  if (neByte %% 4 == 2) {
    x <- x[1:(ndByte - 2)]
  }
  if (neByte %% 4 == 3) {
    x <- x[1:(ndByte - 1)]
  }
  r <- as.raw(x)
  TypeList <- c("logical", "integer", "double", "complex", "character", "raw", "numeric", "int")
  if (!is.character(what) || length(what) != 1 || !(what %in% TypeList)) {
    what <- typeof(what)
  }
  if (what == "raw") {
    return(r)
  }
  if (is.na(size)) {
    size <- switch(
      match(what, TypeList),
      4,
      4,
      8,
      16,
      2,
      1,
      8,
      4
    )
  }
  if (what == "character") {
    rlen <- size * ceiling(length(r) / size)
    length(r) <- rlen
  }
  n <- length(r)
  if (n %% size) {
    stop("raw2bin: number of elements in 'r' is not multiple of 'size'")
  }
  x <- readBin(r, what, n = n %/% size, size = size, signed = signed, endian = endian)
  if (what == "character") {
    x <- paste(x, collapse = "")
  }
  return(x)
}

#' @rdname base64decode
#' @param x The value vector to be encoded.
#' @export
#' @examples
#' \dontrun{
#' # you need to have the bitops package installed to run the example
#' x <- c(10, 0.2, 123456)
#' (z <- base64encode(x = x))
#' base64decode(z = z, what = "numeric")
#' (x <- as.integer(x))
#' (z <- base64encode(x = x))
#' base64decode(z = z, what = "int")
#' }
base64encode <- function(x, size = NA, endian = .Platform$endian) {
  
  verify_suggested(c("bitops"))
  
  if ((typeof(x) != "character") & (typeof(x) != "raw")) {
    x <- writeBin(x, raw(), size = size, endian = endian)
  }
  if ((typeof(x) == "character") & (typeof(x) != "raw")) {
    nlen <- nchar(x)
    x <- writeBin(x, raw(), size = size, endian = endian)
    length(x) <- nlen
  }
  x <- as.integer(x)
  ndByte <- length(x)
  nBlock <- ceiling(ndByte / 3)
  neByte <- 4 * nBlock
  if (ndByte < 3 * nBlock) {
    x[(ndByte + 1):(3 * nBlock)] <- 0
  }
  dim(x) <- c(3, nBlock)
  y <- matrix(as.integer(0), 4, nBlock)
  y[1, ] <- bitops::bitShiftR(x[1, ], 2)
  y[2, ] <- bitops::bitOr(bitops::bitShiftL(x[1, ], 4), bitops::bitShiftR(x[2, ], 4))
  y[3, ] <- bitops::bitOr(bitops::bitShiftL(x[2, ], 2), bitops::bitShiftR(x[3, ], 6))
  y[4, ] <- x[3, ]
  y <- bitops::bitAnd(y, 63)
  alpha <- "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"
  alpha <- strsplit(alpha, NULL)[[1]]
  z <- alpha[y + 1]
  npbytes <- 3 * nBlock - ndByte
  if (npbytes > 0) {
    z[(neByte - npbytes + 1):neByte] <- "="
  }
  z <- paste(z, collapse = "")
  return(z)
}