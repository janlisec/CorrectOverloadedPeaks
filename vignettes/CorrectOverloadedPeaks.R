## ----preload------------------------------------------------------------------
library(CorrectOverloadedPeaks)
data("mzXML_data")

## -----------------------------------------------------------------------------
pk <- CorrectOverloadedPeaks::ModelGaussPeak(height=10^7, width=3, scan_rate=10, e=0, ds=8*10^6, base_line=10^2)
plot(pk, main="Gaussian peak of true intensity 10^7 but cutt off at 8*10^6")

## -----------------------------------------------------------------------------
idx <- pk[,"int"]>0.005 * max(pk[,"int"])
tmp <- CorrectOverloadedPeaks::FitGaussPeak(x=pk[idx,"rt"], y=pk[idx,"int"], silent=FALSE, xlab="RT", ylab="Intensity")

## -----------------------------------------------------------------------------
tmp <- CorrectOverloadedPeaks::CorrectOverloadedPeaks(data=mzXML_data, method="EMG", testing=TRUE)

## -----------------------------------------------------------------------------
load("cor_df_all.RData")
head(cor_df_all[[1]][[1]])
tmp <- CorrectOverloadedPeaks::FitPeakByIsotopicRatio(cor_df=cor_df_all[[1]][[1]], silent=FALSE)

## -----------------------------------------------------------------------------
tmp <- CorrectOverloadedPeaks::FitGaussPeak(x=cor_df_all[[1]][[1]][,"RT"], y=cor_df_all[[1]][[1]][,"int0"], silent=FALSE, xlab="RT", ylab="Intensity")

## ---- echo=FALSE, prompt=FALSE------------------------------------------------
if(file.exists("cor_df_all.RData")) file.remove("cor_df_all.RData")
if(file.exists("S5_35_01_2241_Int+LM.mzXML.pdf")) file.remove("S5_35_01_2241_Int+LM.mzXML.pdf")
if(file.exists("mzXML_data.pdf")) file.remove("mzXML_data.pdf")

