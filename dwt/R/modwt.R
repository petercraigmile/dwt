## ======================================================================
## Copyright 1998--2008, Peter F. Craigmile, All Rights Reserved
## Address comments about this software to pfc@stat.osu.edu.
##
## GNU GENERAL PUBLIC LICENSE, Version 3
## https://www.gnu.org/licenses/gpl-3.0.txt
## ======================================================================


## ======================================================================
## File     : modwt.R
## Contains : Defining the maximum overlap discrete wavelet
##            transform (MODWT) and its inverse (IMODWT) in R.
## Version  : 1.0
## Updated  : pfc@stat.osu.edu, Jul 2007.
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




modwt <- function (x, wavelet = "LA8",
                   nlevels = floor(log2(length(x)))) {
  ## ======================================================================
  ## Purpose : Perform a MODWT decomposition of the data 'x' to
  ##           level 'nlevels' using a given 'wavelet'.
  ## ======================================================================
  
  filter <- dwt.filter(wavelet)
  
  js <- 1:nlevels
  
  N <- length(x)
  Nj <- N/2^js
  L  <- length(filter$scaling)
  Lj       <- (2^js-1)*(filter$L-1)+1
  Bj       <- Lj-1
  Bj[Bj>N] <- N
  
  V <- x
  W <- list(nlevels)
  
  for (j in js) {
    results <- modwt.cascade.next(V, j, filter)
    V      <- results$V
    W[[j]] <- results$W
  }
  
  structure(list(x = x, W = W, V = V, nlevels = nlevels, N = N,
                 Lj=Lj, Mj=N-Lj+1, Bj=Bj,
                 filter = filter),
            class = "modwt")
}




imodwt <- function (dx, use.C=TRUE) {
  ## ======================================================================
  ## Purpose : Performs an inverse MODWT of the modwt class 'dx'.
  ## ======================================================================
  
  V <- dx$V
  for (j in dx$nlevels:1)
    V <- modwt.cascade.previous(V, dx$W[[j]], j, dx$filter, use.C)
  
  V
}



modwt.wavelets.zero <- function (dx, levels) {
  ## ======================================================================
  ## Purpose : Make the wavelet coefficients on level 'levels' of the
  ##           modwt class 'dx' zero.
  ## ======================================================================

  dz <- dx
  for (j in levels)
    dz$W[[j]] <- rep(0,length(dz$W[[j]]))
  dz
}




## ======================================================================
## Purpose: Calculates the level 'j' detail coefficients of the MODWT
##          object 'dx'.
## Assumes: 1 <= j <= dx$nlevels.
## History: pfc@stat.osu.edu, Jul 2002.
## ======================================================================

modwt.detail <- function (dx, j, use.C=TRUE)
{
  dz   <- modwt.wavelets.zero(dx, seq(dx$nlevels)[-j])
  dz$V <- dz$V * 0
  imodwt(dz, use.C)
}




## ======================================================================
## Purpose: Calculates the level 'j' smooth coefficients of the MODWT
##          object 'dx'.
## Assumes: 1 <= j <= dx$nlevels.
## History: pfc@stat.osu.edu, Jul 2002.
## ======================================================================

modwt.smooth <- function (dx, j, use.C=TRUE)
{
  imodwt(modwt.wavelets.zero(dx, 1:j), use.C)
}


