## DWT

Email questions about the code to pfc &lt;AT&gt; stat.osu.edu

R routines for carrying out the Discrete Wavelet Transforms (DWT) and Maximum Overlap Discrete Wavelet Transforms (MODWT) for univariate time series.

To install the R package make sure that you have the `devtools` R package installed first and then type
```
devtools::install_github("dwt", user="petercraigmile", subdir="dwt") 
```

The R code contains functions to carry out the wavelet transforms in R and in R with C (to speed up the calculations).  I need to check the two versions of the MODWT code.

The R package is contained in the `dwt` folder.

The documentation for this R package needs more work.



### References:

D. B. Percival and A. T. Walden. Wavelet Methods for Time Series Analysis. Cambridge University
Press, Cambridge, 2000.