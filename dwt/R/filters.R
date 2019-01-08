## ======================================================================
## Copyright 1998--, Peter F. Craigmile, All Rights Reserved
## Address comments about this software to pfc@stat.osu.edu.
##
## GNU GENERAL PUBLIC LICENSE, Version 3
## https://www.gnu.org/licenses/gpl-3.0.txt
## ======================================================================


## ======================================================================
## File     : filters.R
## Contains : Some standard operations for filters.
## Version  : 1.1
## Updated  : pfc@stat.osu.edu, Apr 2009.
##
## References:
##
##     D. B. Percival and A. T. Walden (2000),
##     Wavelet Methods for Time Series Analysis.
##     Cambridge, England: Cambridge University Press.
## ======================================================================




downsample.by.two <- function (x, even=TRUE) {
  ## ======================================================================
  ## Return the even indexes of 'x' (if 'even' is TRUE) or the odd indexes
  ## of 'x' (if 'even' is FALSE)  
  ## ======================================================================

  if (even)
    x[seq(2, length(x), 2)]
  else
    x[seq(1, length(x), 2)]
}



upsample.by.two <- function (x, before=FALSE) {
  ## ======================================================================
  ## Interleave one zero between the values of 'x'.
  ## If before is TRUE, each value of 'x' occurs befores the zero.
  ## ======================================================================

  if (before)
    as.numeric(rbind(x, 0))
  else
    as.numeric(rbind(0, x))
}



upsample.by <- function (x, n, before=TRUE) {
  ## ======================================================================
  ## Interleave 'n' zeroes between each values of 'x'.
  ## If before is TRUE, each value of 'x' occurs befores the zeroes
  ## ======================================================================

    if (before) {
        
        as.numeric(rbind(x, matrix(0, n, length(x))))
        
    } else {
        
        as.numeric(rbind(matrix(0, n, length(x)), x))
    }
}



periodize.filter <- function (f, n, L=length(f)) {
  ## ======================================================================
  ## periodize the filter 'f' of length 'L'
  ## into a vector of length 'n'.
  ## ======================================================================

  ## This method is slower
  ## as.numeric(rowsum(f, (0:(L-1)) %% n))

  ## use this instead
  rowSums(matrix(c(f, rep(0, n - L %% n)), n))
}




circularly.filter <- function (x, f, L=length(f)) {
  ## ======================================================================
  ## Circularly filter the sequence 'x' using the filter 'f'.
  ## ======================================================================
  
  n <- length(x)
  
  ## Do we need to periodize the filter?
  if (L > n)
    f <- periodize.filter(f, n, L)

  filter(x, f, "convolution", sides=1, circular=TRUE)
}



center.of.energy <- function (f) {
  ## ======================================================================
  ## Calculate the center of energy of the filter 'f'
  ## References: Wickerhauser (1994, p.171 and p.341)
  ##             Percival and Walden (2000, p.118)
  ## ======================================================================

  f.sq <- f^2
  sum((0:(length(f)-1)) * f.sq) / sum(f.sq)
}
