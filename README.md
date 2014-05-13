## DWT

Email questions about the code to pfc &lt;AT&gt; stat.osu.edu

R routines for carrying out the Discrete Wavelet Transforms (DWT) and Maximum Overlap Discrete Wavelet Transforms (MODWT) for univariate time series.

To install the R package make sure that you have the `devtools` R package installed first and then type
```
devtools::install_github("dwt", user="petercraigmile", subdir="dwt") 
```

The R code contains functions to carry out the wavelet transforms in R and in R with C (to speed up the calculations).


The R package is contained in the `dwt` folder.




### Installation

For Mac OS and Linux, make sure you have the C and fortran compilers
installed.  After installing the `devtools` R package type

```
devtools::install_github("dwt", user="petercraigmile", subdir="dwt") 
```

For windows, download the file <a href="https://github.com/petercraigmile/dwt/raw/master/releases/current/dwt.zip">dwt.zip</a> from the `releases` folder.    Then:

1. Open the R gui, by doubling clicking on the R icon.

2. Click on the `Packages` menu and select `Install package(s) from local zip files...`.  Find the zip file and press `Open`.

3. Your R package will be installed.


After installation, type `library(dwt)` to use the R library.

### Issues

1. I need to check the two versions of the MODWT code.

2. The documentation for this R package needs more work.


### References:

D. B. Percival and A. T. Walden. Wavelet Methods for Time Series Analysis. Cambridge University
Press, Cambridge, 2000.