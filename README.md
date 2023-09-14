
# CorrectOverloadedPeaks

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/CorrectOverloadedPeaks)](https://CRAN.R-project.org/package=CorrectOverloadedPeaks)
<!-- badges: end -->

Time series data are often analysed for peak signals. Mass spectrometry data may contain 
flat top peaks due to technical limitations (i.e. detector saturation, DS). Flat top peaks 
can also be termed 'overloaded' signals. Extracting the peak height to infer signal 
intensity will obviously give wrong results for flat top peaks. However, using the peak 
shape in the non-distorted fraction of the signal (intensity below DS), the true peak
shape can be modeled mathematically. 
This modelling is the core task of `CorrectOverloadedPeaks`. The R package accepts data in 
xcmsRaw and mzXML format as input. Overloaded signals are detected automatically and modified 
using an Gaussian or Isotopic-Ratio approach, QC plots are generated and corrected data are 
stored within the original xcmsRaw or mzXML respectively to allow further processing.
This way `CorrectOverloadedPeaks` can be incorporated in any metabolomics pipeline. Some
utility functions are additionally exported, i.e. `read.mzData()` and `FitGaussPeak()`.

## Installation

You can install the development version of CorrectOverloadedPeaks from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("janlisec/CorrectOverloadedPeaks")
```

or install from CRAN otherwise.

## Example

This is a basic example, modelling a flat topped peak first and restoring the true shape 
assuming a Gaussian peak shape afterwards.

``` r
pk <- CorrectOverloadedPeaks::ModelGaussPeak(height=10^7, width=3, scan_rate=10, e=0, ds=8*10^6, base_line=10^2)
plot(pk, main="Gaussian peak of true intensity 10^7 but cutt off at 8*10^6")
idx <- pk[,"int"]>0.005 * max(pk[,"int"])
tmp <- CorrectOverloadedPeaks::FitGaussPeak(x=pk[idx,"rt"], y=pk[idx,"int"], silent=FALSE, xlab="RT", ylab="Intensity")
```
Next, we load some real life measurement data and correct the two overloaded peaks contained.

``` r
data("mzXML_data", package = "CorrectOverloadedPeaks")
tmp <- CorrectOverloadedPeaks::CorrectOverloadedPeaks(data=mzXML_data, method="EMG", testing=TRUE)
```
