## DWT

Email questions about the code to pfc &lt;AT&gt; stat.osu.edu

R routines for carrying out the Discrete Wavelet Transforms (DWT) and
Maximum Overlap Discrete Wavelet Transforms (MODWT) for univariate
time series.   C code versions are included to speed up calculations,
and are used by defauly.

The R package is contained in the `dwt` folder.




### Installation

Please make sure you have a C compiler installed for your version of R. Then, after installing the devtools R package, you can install this R package using

```
devtools::install_github("petercraigmile/dwt/dwt") 
```

After installation, type `library(dwt)` to use the R library.



### Issues

1. The plot command is not working.

2. The documentation for this R package needs more work.


### References:

D. B. Percival and A. T. Walden. Wavelet Methods for Time Series Analysis. Cambridge University
Press, Cambridge, 2000.
