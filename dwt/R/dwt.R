## ======================================================================
## Copyright 1998--2009, Peter F. Craigmile, All Rights Reserved
## Address comments about this software to pfc@stat.osu.edu.
##
## GNU GENERAL PUBLIC LICENSE, Version 3
## https://www.gnu.org/licenses/gpl-3.0.txt
## ======================================================================


## ======================================================================
## File     : dwt.R
## Contains : Defining the discrete wavelet transform (DWT) and its
##            inverse (IDWT) in R.
## Version  : 1.01
## Updated  : pfc@stat.osu.edu, Apr 2008.
## Requires : 'dwt_filters.R' and 'filters.R'
##
## Notes    : This R code does not require extra dynamically allocated
##            C code, and so is very portable, but several times slower
##            than dedicated C code would be.
##            I could have optimized the code more, but I left the
##            filtering details intact for educational purposes.
##
## References:
##
##     D. B. Percival and A. T. Walden (2000),
##     Wavelet Methods for Time Series Analysis.
##     Cambridge, England: Cambridge University Press.
## ======================================================================



dwt <- function (x, wavelet = "LA8",
                 nlevels = floor(log2(length(x))), use.C=TRUE) {
  ## ======================================================================
  ## Purpose : Perform a DWT decomposition of the data 'x' to
  ##           level 'nlevels' using a given 'wavelet'.
  ## Assumes : Can compose to 'nlevels' levels;
  ##           i.e., length of 'x' is divisible by 2^'nlevels'. 
  ## ======================================================================
  
  filter <- dwt.filter(wavelet)
  
  js <- 1:nlevels
  
  N <- length(x)
  Nj <- N/2^js
  L  <- length(filter$scaling)
  Bj <- ceiling((L - 2) * (1 - 2^(-js)))
  Lj <- (2^js - 1) * (L - 1) + 1
  
  V <- x
  W <- list(nlevels)
  
  for (j in js) {
    results <- dwt.cascade.next(V, filter, use.C)
    V      <- results$V
    W[[j]] <- results$W
  }
  
  structure(list(x = x, W = W, V = V, nlevels = nlevels, N = N, Nj = Nj, 
                 Bj = Bj, Lj = Lj, Mj = Nj - Bj, M = N - sum(Bj),
                 filter = filter),
            class = "dwt")
}




idwt <- function (dx, use.C=TRUE) {
  ## ======================================================================
  ## Purpose : Performs an inverse DWT of the dwt class 'dx'.
  ## ======================================================================
  
  V <- dx$V
  for (j in dx$nlevels:1)
    V <- dwt.cascade.previous(V, dx$W[[j]], dx$filter, use.C)
  
  V
}




dwt.wavelets.zero <- function (dx, levels, ignore.boundary=FALSE) {
  ## ======================================================================
  ## Purpose : Make the wavelet coefficients on level 'levels' of the
  ##           dwt class 'dx' zero.
  ##           If ignore.boundary is TRUE, do not zero out the boundary
  ##           coefficients
  ## ======================================================================

  dz <- dx
  if (ignore.boundary) {
    for (j in levels) {
      if (dx$Mj[j]>0) {
        dz$W[[j]][(dx$Bj[j]+1):dx$Nj[j]] <- rep(0,length(dx$Mj[j]))
      }
    }
  } else {
    for (j in levels)
      dz$W[[j]] <- rep(0,length(dz$W[[j]]))
  }
  dz
}



## ======================================================================
## Purpose: Calculates the level 'j' detail coefficients of the DWT
##          object 'dx'.
## Assumes: 1 <= j <= dx$nlevels.
## History: pfc@stat.osu.edu, Jul 2002.
## ======================================================================

dwt.detail <- function (dx, j, use.C=TRUE)
{
  dz   <- dwt.wavelets.zero(dx, seq(dx$nlevels)[-j])
  dz$V <- dz$V * 0
  idwt(dz, use.C)
}




## ======================================================================
## Purpose: Calculates the level 'j' smooth coefficients of the DWT
##          object 'dx'.
## Assumes: 1 <= j <= dx$nlevels.
## History: pfc@stat.osu.edu, Jul 2002.
## ======================================================================

dwt.smooth <- function (dx, j, use.C=TRUE)
{
  idwt(dwt.wavelets.zero(dx, 1:j), use.C)
}





