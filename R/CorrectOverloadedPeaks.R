#' @title Correct Overloaded Peaks from GC-MS data.
#'
#' @description \code{CorrectOverloadedPeaks} will take an xcmsRaw data structure 
#'     (or any imported mzXML) and search for overloaded peaks within the mass traces. 
#'     It will correct overloaded peaks automatically using an Gaussian or 
#'     IsotopicRatio approach, generate QC plots and write the corrected data 
#'     back into the original xcmsRaw.
#'
#' @details This is a high level function to batch pre-process metabolomics data 
#'     which are partially overloaded before continuing with the standard workflow 
#'     of peak identification etc.. It relies internally on \code{\link{FitGaussPeak}} 
#'     and \code{\link{FitPeakByIsotopicRatio}} to modify data of individual intensity 
#'     signal. Basically the function aims to identify automatically overloaded regions
#'     and extracts base peak chromatograms for all overloaded m/z traces within these 
#'     regions, which are corrected and put back into the original data structure. 
#'     For simplicity some potentially interesting parameters are hidden at the top of 
#'     the function definition. They have been set to values determined empirically to 
#'     be working for a Bruker impact II MS (high-res QTOF) coupled to GC and LC via 
#'     APCI and ESI respectively. For more details please 
#'     see \doi{10.1021/acs.analchem.6b02515}.
#'
#' @references \doi{10.1021/acs.analchem.6b02515}
#'
#' @param data An xcmsRaw-object or an mzXML-object as imported by \code{\link{read.mzXML}}.
#' @param method Either Gauss or EMG (usually better results) or Isoratio (more robust for non-Gaussian peak shapes).
#' @param detection_limit If=1 only peaks hitting detector saturation (ds) will be corrected, can be lowered to 0.95 to catch also peaks going into saturation.
#' @param ds Detector saturation. Will be determined based on data if not specified explicitly.
#' @param silent QC-plots will be generated if silent=FALSE and additional Warnings() will be generated.
#' @param testing Will automatically set silent=FALSE and store all extracted regions with overloaded peaks in the working directory as \code{cor_df_all.RData}.
#' @param attotwm All-the-Time-of-the-World-Mode. If calculation time doesn't matter try this out. :)
#' @param region From an initial QC-Plot file you may reprocess a specific overloaded region. Don't forget to specify the ds parameter explicitly.
#' @param peak You may further restrict the reprocessing to a specific peak within the region.
#'
#' @return
#' An corrected xcmsRaw- or mzXML-object which can be exported to file. Additionally a QC-plot pdf-file if silent=FALSE.
#'
#' @examples
#' \dontrun{
#'   # load mzXML test data
#'   data(mzXML_data)
#'   CorrectOverloadedPeaks(data = mzXML_data, method = "EMG", silent = FALSE)
#' }
#'
#' @seealso \code{\link{ModelGaussPeak}}
#' @seealso \code{\link{FitGaussPeak}}
#' @seealso \code{\link{FitPeakByIsotopicRatio}}
#' @seealso \code{\link{read.mzXML}}
#' @seealso \code{\link{write.mzXML}}
#'
#' @export
#'
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @importFrom graphics legend
#' @importFrom graphics par
#' @importFrom stats median
#'
CorrectOverloadedPeaks <- function(
  data = NULL, 
  method = c("Isoratio", "Gauss", "EMG"), 
  detection_limit = 1, 
  ds = NULL, 
  silent = TRUE, 
  testing = FALSE, 
  attotwm = FALSE, 
  region = NULL, 
  peak = NULL) 
{
  # POTENTIAL PARAMETERS that could be allowed for the user to modify
  # mz_gap_allowed : what is the MS resolution? Mass accuracy in case of overloading? And, more important, expected mass error in case of overloading (0.13==130mDa)
  mz_gap_allowed <- 0.015
  # allowed_scans_missing : get number and position of overloaded regions as well as overloaded masses
  allowed_scans_missing <- 4
  # mz_dev : to test for mass stability [in ppm]
  mz_dev <- 15
  # corr_scan_offset : what is the upstream/downstream offset used for mass/int correction
  # depending on the scan rate and the separation quality/slope of peak front/tail
  corr_scan_offset <- c(-10, 10)
  # ds_cutoff : lower level for iso-ratio determination; upper level for non-linearity/correction-assumption
  ds_cutoff <- c(0.005, 0.9)
  # cor_fac_front : for 'isoratio': use only front (TRUE) or tail (FALSE) information to correct data
  cor_fac_front <- TRUE
  # limit_of_correction [time unit of sample]: peaks which are overloaded for more than this should not attempted to be corrected
  limit_of_correction <- 8
  # allowed_front_tail_ratio : what is the maximum ratio allowed in min(front)/min(tail); if aftr=10 than fronting is assumed for 10 fold difference; set to Inf to avaoid peak shape testing
  aftr <- 2.5
  # try to correct for ion suppression
  correct_ion_suppression <- FALSE
  # in Gauss method cut off baseline to improve fitting
  account_for_baseline_offset <- TRUE

  if (testing) silent <- FALSE

  method <- match.arg(method)
  
  # internal functions
  DetermineDetectorSaturation <- function(ds, pks, inf) {
    if (is.null(ds)) {
      ds <- max(sapply(pks, function(x) {
        max(x[, 2])
      }))
      c1 <- sum(sapply(pks, function(x) {
        sum(x[, 2] == ds)
      })) >= 3
      c2 <- sum(table(unlist(sapply(pks, function(x) {
        x[x[, 2] > (ds / 2), 2]
      }))) >= sum(sapply(pks, function(x) {
        sum(x[, 2] == ds)
      })) / 2) == 1
      c3 <- all(sapply(pks, function(x) {
        max(x[, 2]) == ds
      }) == (inf[, "basePeakIntensity"] == ds))
      overloaded_peaks_present <- c1 && c2 && c3
    } else {
      overloaded_peaks_present <- any(sapply(pks, function(x) {
        any(x[, 2] >= ds)
      }))
    }
    attr(ds, "opp") <- overloaded_peaks_present
    return(ds)
  }
  GetOverloadedRegions <- function(inf, pks, lim, allowed_scans_missing, mz_gap_allowed) {
    tmp <- GroupByGaps(x = inf[which(inf[, "basePeakIntensity"] >= lim), "retentionTime"], gap = allowed_scans_missing * median(diff(inf[, "retentionTime"])))
    tmp <- split(inf[which(inf[, "basePeakIntensity"] >= lim), c("seqNum", "retentionTime", "basePeakMZ")], tmp) # ; length(ovreg)
    mzaff_all <- lapply(1:length(tmp), function(i) {
      sapply(tmp[[i]][, "seqNum"], function(j) {
        pks[[j]][which(pks[[j]][, 2] >= lim), 1]
      })
    })
    mzovld <- lapply(mzaff_all, function(x) {
      sapply(split(unlist(x), GroupByGaps(unlist(x), mz_gap_allowed)), median)
    })
    for (i in 1:length(mzovld)) attr(mzovld[[i]], "seqNum") <- range(tmp[[i]][, "seqNum"])
    return(mzovld)
  }
  ExtractRawData <- function(mz, scans, pks, lim, mz_gap_allowed) {
    # extract intensity of M+0 and M+x until the first isotope no exceeding ds in wide window
    idx_iso <- 1
    cor_df <- data.frame(
      inf[scans, c("seqNum", "retentionTime")],
      t(sapply(scans, function(j) {
        pks[[j]][which.min(abs(mz - pks[[j]][, 1])), , drop = F]
      })),
      t(sapply(scans, function(j) {
        pks[[j]][which.min(abs(mz + idx_iso * 1.0033 - pks[[j]][, 1])), , drop = F]
      }))
    )
    colnames(cor_df) <- c("Scan", "RT", "mz0", "int0", "mz1", "int1")
    check_iso_int <- any(cor_df[, "int1"] >= detection_limit * ds)
    while (check_iso_int) {
      idx_iso <- idx_iso + 1
      cor_df <- cbind(cor_df, data.frame(t(sapply(scans, function(j) {
        pks[[j]][which.min(abs(mz + idx_iso * 1.0033 - pks[[j]][, 1])), , drop = F]
      }))))
      colnames(cor_df)[ncol(cor_df) - c(1, 0)] <- paste0(c("mz", "int"), idx_iso)
      check_iso_int <- any(cor_df[, paste0("int", idx_iso)] >= detection_limit * ds)
    }
    cor_df <- cbind(cor_df, "modified" = FALSE)
    # restrict to overloaded area for this specific M0
    idx <- range(which(cor_df[, "int0"] >= (ds_cutoff[2] * detection_limit * ds)))
    if ((idx[1] + corr_scan_offset[1]) >= 1) {
      idx[1] <- idx[1] + corr_scan_offset[1]
    } else {
      idx[1] <- 1
    }
    if ((idx[2] + corr_scan_offset[2]) <= nrow(cor_df)) {
      idx[2] <- idx[2] + corr_scan_offset[2]
    } else {
      idx[2] <- nrow(cor_df)
    }
    idx <- seq(idx[1], idx[2])
    cor_df <- cor_df[idx, ]
    # test if this area contains two peaks (e.g. isoforms) of same mz and limit to larger overloaded and give a warning
    # tmp <- cor_df[,"int0"]>(ds_cutoff[2]*detection_limit*ds)
    # diff(cumsum(tmp))
    # remove scans with strong mz deviation for M0
    f <- abs(cor_df[, "mz0"] - mz) < mz_gap_allowed
    if (any(!f)) cor_df <- cor_df[cumsum(f) > 0, ]
    # limit to interesting area
    idx <- range(which(cor_df[, "int0"] > lim))
    cor_df <- cor_df[seq(idx[1], idx[2]), ]
    return(cor_df)
  }
  PutDataBack <- function(pks, cor_df) {
    if (any(cor_df[, "modified"])) {
      for (j in cor_df$Scan[cor_df[, "modified"]]) {
        idx <- which(cor_df$Scan == j)
        pks[[j]][which.min(abs(pks[[j]][, 1] - cor_df[idx, "mz0"])), 2] <- cor_df[idx, "int0"]
        pks[[j]][which.min(abs(pks[[j]][, 1] - cor_df[idx, "mz1"])), 2] <- cor_df[idx, "int1"]
        pks[[j]][which.min(abs(pks[[j]][, 1] - cor_df[idx, "mz2"])), 2] <- cor_df[idx, "int2"]
      }
    }
    return(pks)
  }

  # extract peak and scan infos
  if (inherits(data, "xcmsRaw")) {
    scannum <- as.integer(diff(c(data@scanindex, length(data@env$mz))))
    pks <- sapply(split(matrix(c(data@env$mz, data@env$intensity), ncol = 2), factor(rep(1:length(scannum), times = scannum))), matrix, ncol = 2)
    # $$JL: fixed in V1.2.17$$
    # inf <- data.frame("seqNum"=data@acquisitionNum,"retentionTime"=data@scantime,"basePeakMZ"=sapply(pks,function(x){x[which.max(x[,2])[1],1]}),"basePeakIntensity"=sapply(pks,function(x){max(x[,2])}))
    inf <- data.frame("seqNum" = 1:length(data@scanindex), "retentionTime" = data@scantime, "basePeakMZ" = sapply(pks, function(x) {
      x[which.max(x[, 2])[1], 1]
    }), "basePeakIntensity" = sapply(pks, function(x) {
      max(x[, 2])
    }))
  }

  if (inherits(data, "mzXML")) {
    pks <- lapply(data[["scan"]], function(x) {
      as.matrix(data.frame("mz" = as.numeric(x$mass), "intensity" = as.numeric(x$peaks)))
    })
    tmp.rt <- sapply(data[["scan"]], function(x) {
      x <- strsplit(x$scanAttr, " ")[[1]]
      x[grep("retentionTime", x)]
    })
    tmp.rt <- as.numeric(gsub("S\"", "", gsub("retentionTime=\"PT", "", tmp.rt)))
    inf <- data.frame("seqNum" = 1:length(pks), "retentionTime" = tmp.rt, "basePeakMZ" = sapply(pks, function(x) {
      x[which.max(x[, 2])[1], 1]
    }), "basePeakIntensity" = sapply(pks, function(x) {
      max(x[, 2])
    }))
  }

  # if an attribute 'file_in' is attached to 'data' it will be used for warnings and plot names
  if (!is.null(attr(data, "file_in"))) file_in <- attr(data, "file_in") else file_in <- ""

  # determine value of detector saturation (if not specified explicitly) and presence of overloaded peaks
  ds <- DetermineDetectorSaturation(ds = ds, pks = pks, inf = inf)

  if (!silent) cat(paste("\nProcessing...", file_in, "\n"))

  # apply correction if overloaded peaks are present
  if (!attr(ds, "opp")) {
    if (!silent) print("Found no overloaded peaks in sample. Return data as is.")
  } else {
    # file name for plot output
    if (!silent) grDevices::pdf(ifelse(file_in == "", paste0(format(Sys.time(), "%Y-%m-%d_%H%M%S"), "_CorrectOverloadedPeaks_ControlPlot.pdf"), paste0(file_in, ".pdf")))

    # search rt-regions and masses above detector saturation
    mzovld <- GetOverloadedRegions(inf = inf, pks = pks, lim = detection_limit * ds, allowed_scans_missing = allowed_scans_missing, mz_gap_allowed = mz_gap_allowed)
    if (!silent) cat(paste("\nTrying to correct", length(mzovld), "overloaded regions.\n"))

    # write all extracted raw data into list for individual file for testing purposes
    if (testing) {
      cor_df_all <- lapply(1:length(mzovld), function(i) {
        vector("list", length = length(mzovld[[i]]))
      })
    }

    # correct overloaded intensity data in 'i' regions by extrapolating based on isotopic ratio or gaussian approximation
    for (i in 1:length(mzovld)) {
      tmp_rng <- attr(mzovld[[i]], "seqNum")
      if ((tmp_rng[1] + corr_scan_offset[1]) >= 1) tmp_rng[1] <- tmp_rng[1] + corr_scan_offset[1]
      if ((tmp_rng[2] + corr_scan_offset[2]) <= nrow(inf)) tmp_rng[2] <- tmp_rng[2] + corr_scan_offset[2]
      scans <- seq(min(tmp_rng), max(tmp_rng))

      # and all M0 within these regions
      for (k in 1:length(mzovld[[i]])) {
        # extract data/put M0-M2 infos together
        cor_df <- ExtractRawData(mz = mzovld[[i]][k], scans = scans, pks = pks, lim = ds_cutoff[1] * ds, mz_gap_allowed = mz_gap_allowed)

        # keep next line for testing
        if (testing) cor_df_all[[i]][[k]] <- cor_df

        # error checks/warnings
        # [$$] is currently not used (to many warnings)
        if (!silent & FALSE) {
          # test for mass stability
          if (!all(apply(cor_df[, c("mz0", "mz1", "mz2")], 2, function(x) {
            diff(range(x)) < mz_dev * median(x) / 10^6
          }))) {
            warning(paste(basename(file_in), "region", i, "m/z", k, "Strong mass diff"))
          }
          # test if M1 and M2 are overloaded
          if (all(apply(cor_df[, c("int1", "int2")], 2, function(x) {
            max(x) > ds_cutoff[2] * ds
          }))) {
            warning(paste(basename(file_in), "region", i, "m/z", k, "Both isotopes approach detector saturation"))
          }
          # test if M1 and M2 are flatened
          if (all(apply(cor_df[, c("int1", "int2")], 2, is.FlatTopPeak))) {
            warning(paste(basename(file_in), "region", i, "m/z", k, "Both isotopes are flat-toped (while below detector saturation)"))
          }
        }

        # get the scans which need to be corrected, this is usually determined by 'ds_cutoff' and 'detection_limit'
        idx <- range(which(cor_df[, "int0"] >= (ds_cutoff[2] * detection_limit * ds)))
        idx <- seq(idx[1], idx[2])

        # however, in regions with several overloaded signals we observe ion suppression for 'smaller' peaks which could be isotopes, fragments
        # but also M+H in case that a fragment is the base peak. Therefore, we may want to modify the cut_off parameter in FitGaussPeak to not include
        if (correct_ion_suppression) {
          # data points already affected by ion suppression. We can test that by:
          test_ionsupp <- attr(mzovld[[i]], "seqNum") - range(cor_df[, "Scan"])
          test_ionsupp <- test_ionsupp[1] < 4 | test_ionsupp[2] > -4
          # test_ionsupp <- test_ionsupp[1]<4 & test_ionsupp[2]>-4
          if (test_ionsupp) {
            # and adjust idx accordingly
            cutoff_modifier <- 0.4
            idx <- range(which(cor_df[, "int0"] >= (cutoff_modifier * ds_cutoff[2] * detection_limit * ds)))
            idx <- seq(idx[1], idx[2])
          }
        } else {
          test_ionsupp <- FALSE
        }

        # check for various potential errors and drawbacks
        err_code <- 0
        # set a limit to meaningfull correction
        if (diff(range(cor_df[idx, "RT"])) >= limit_of_correction) err_code <- 1
        # check if front/tail is missing (may happen, if import window was to narrow)
        if (min(idx) == 1) err_code <- 2
        if (max(idx) == nrow(cor_df)) err_code <- 3

        # temporary error code for testing
        # specifying i and k as identified from a QC-plot file here allows to interrupt at this specific peak
        if (!is.null(region)) {
          err_code <- ifelse(i != region, 9, err_code)
        }
        if (!is.null(peak)) {
          err_code <- ifelse(k != peak, 9, err_code)
        }

        # print some error messages
        if (err_code > 0) {
          tmp.rt <- round(mean(cor_df[idx, "RT"]), 2)
          tmp.mz <- round(mean(cor_df[idx, "mz0"]), 4)
          if (err_code == 1) print(paste("Did not correct a peak with overloading wider than", limit_of_correction, "at RT =", tmp.rt, "for mz =", tmp.mz))
          if (err_code == 2) print(paste("Did not correct a peak without detectable front at RT =", tmp.rt, "for mz =", tmp.mz))
          if (err_code == 3) print(paste("Did not correct a peak without detectable tail at RT =", tmp.rt, "for mz =", tmp.mz))
          if (err_code == 9) print("Skip peak because 'region' and/or 'peak' are specified.")
        } else {
          # browser()
          if (!silent) print(paste("Processing Region/Mass:", paste(i, "/", k)))

          # calculate IsotopicRatioFit as basis or final, correct only M0 intensity because M+1 will be corrected independently if overloaded
          cor_df_mod <- cor_df
          cor_df_mod[, "int0"] <- FitPeakByIsotopicRatio(cor_df, idx = idx, silent = silent)[, "int0"]

          # calculate gauss fit aditionally if specified
          if (method == "Gauss" | method == "EMG") {
            front_scans <- 1:ifelse(min(idx) <= 1, 1, min(idx) - 1)
            tail_scans <- ifelse(max(idx) < nrow(cor_df), max(idx) + 1, nrow(cor_df)):nrow(cor_df)
            # is peak tail or front skewed??
            # assumption: intensity should start and fall back to baseline of similar vlaue, if not probably fronting/tailing occur
            test_peak_shape <- min(cor_df[front_scans, "int0"]) / min(cor_df[tail_scans, "int0"])
            gauss_scale <- c(0.8, length(idx) * 2)
            gauss_steps <- length(idx) * 10
            gauss_main <- paste("mz =", round(median(cor_df[, "mz0"]), 4), paste0("(Region ", i, ", Peak ", k, ")"))
            # fronting
            if (test_peak_shape > aftr) {
              gauss_fit <- FitGaussPeak(x = cor_df[, "RT"], y = cor_df[, "int0"], scale_range = gauss_scale, steps = gauss_steps, idx = idx, strip_data = "front", account_for_baseline_offset = account_for_baseline_offset, method = method, silent = silent, xlab = "RT", ylab = "Int", main = gauss_main)
              if (!silent) graphics::legend(x = "bottomleft", "fronting detected", bty = "n", text.col = 2)
            }
            # tailing
            if (test_peak_shape < (1 / aftr)) {
              gauss_fit <- FitGaussPeak(x = cor_df[, "RT"], y = cor_df[, "int0"], scale_range = gauss_scale, steps = gauss_steps, idx = idx, strip_data = "tail", account_for_baseline_offset = account_for_baseline_offset, method = method, silent = silent, xlab = "RT", ylab = "Int", main = gauss_main)
              if (!silent) graphics::legend(x = "bottomright", "tailing detected", bty = "n", text.col = 2)
            }
            # more or less symmetric peak
            if (test_peak_shape >= (1 / aftr) & test_peak_shape <= aftr) {
              gauss_fit <- FitGaussPeak(x = cor_df[, "RT"], y = cor_df[, "int0"], scale_range = gauss_scale, steps = gauss_steps, idx = idx, method = method, account_for_baseline_offset = account_for_baseline_offset, silent = silent, xlab = "RT", ylab = "Int", main = gauss_main)
            }
            if (!silent & test_ionsupp) graphics::legend(x = "bottom", "ion suppression detected", bty = "n", text.col = 2)
            # if isoratio-fit is not higher than gauss solution (should not happen) then substitute and keep iso-solution elsewise as a fallback
            if (max(cor_df_mod[, "int0"]) < max(gauss_fit)) cor_df_mod[, "int0"] <- gauss_fit
          }
          if (attotwm) {
            # re-initialize index
            idx <- range(which(cor_df[, "int0"] >= (ds_cutoff[2] * detection_limit * ds)))
            idx <- seq(idx[1], idx[2])
            # generate an overview for results of various parameter combinations
            attotwm_m <- rep(c("Gauss", "EMG"), each = 6)
            attotwm_c <- rep(c(0.4, 0.9, 0.4, 0.9), each = 3)
            attotwm_s <- rep(c("front", "none", "tail"), times = 4)
            attotwm_r <- NULL
            graphics::par(mfrow = c(4, 3))
            graphics::par(cex = 0.4)
            for (attotwm_i in 1:12) {
              attotwm_r <- c(attotwm_r, max(FitGaussPeak(
                x = cor_df[, "RT"], y = cor_df[, "int0"], scale_range = gauss_scale, steps = gauss_steps, account_for_baseline_offset = account_for_baseline_offset, silent = silent, xlab = "RT", ylab = "Int",
                cutoff = attotwm_c[attotwm_i],
                strip_data = attotwm_s[attotwm_i],
                method = attotwm_m[attotwm_i],
                main = paste(attotwm_m[attotwm_i], attotwm_c[attotwm_i], attotwm_s[attotwm_i], sep = ", ")
              )))
            }
            graphics::par(mfrow = c(1, 1))
            graphics::par(cex = 1)
            gauss_fit <- FitGaussPeak(
              x = cor_df[, "RT"], y = cor_df[, "int0"], scale_range = gauss_scale, steps = gauss_steps, account_for_baseline_offset = account_for_baseline_offset, silent = silent, xlab = "RT", ylab = "Int",
              cutoff = attotwm_c[which.max(attotwm_r)],
              strip_data = attotwm_s[which.max(attotwm_r)],
              method = attotwm_m[which.max(attotwm_r)],
              main = paste(attotwm_m[which.max(attotwm_r)], attotwm_c[which.max(attotwm_r)], attotwm_s[which.max(attotwm_r)], sep = ", ")
            )
            if (max(cor_df_mod[, "int0"]) < max(gauss_fit)) cor_df_mod[, "int0"] <- gauss_fit
          }
          cor_df_mod[idx, "modified"] <- TRUE

          # put corrected data back into temporary peak list
          pks <- PutDataBack(pks = pks, cor_df = cor_df_mod)
        }
      }
    }

    # close plotting device
    if (!silent) grDevices::dev.off()

    # assign to base
    if (testing) {
      print("Storing non-corrected data information in 'cor_df_all.RData'")
      save(cor_df_all, file = "cor_df_all.RData")
    }

    # put corrected data back into file
    if (inherits(data, "xcmsRaw")) {
      # data@env$mz <- unlist(sapply(pks,function(x){x[,1]}))
      data@env$intensity <- unlist(sapply(pks, function(x) {
        x[, 2]
      }))
    }
    if (inherits(data, "mzXML")) {
      # for (i in 1:length(pks)) data[["scan"]][[i]][["mass"]]  <- pks[[i]][,1]
      # print(table(sapply(1:length(pks), function(i) {all(data[["scan"]][[i]][["peaks"]]==round(pks[[i]][,2]))})))
      for (i in 1:length(pks)) data[["scan"]][[i]][["peaks"]] <- round(pks[[i]][, 2])
      # print(table(sapply(1:length(pks), function(i) {all(data[["scan"]][[i]][["peaks"]]==round(pks[[i]][,2]))})))
    }
  }
  invisible(data)
}
