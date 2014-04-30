## ======================================================================
## File     : modwt_time_freq.R
## Purpose  : R code for aligning wavelet and scaling coefficients with time
##            or approximately with frequency in the MODWT
## Updated  : pfc@stat.osu.edu, Jun 2009.
## 
## Copyright 2002--, Peter F. Craigmile, All Rights Reserved
## Address comments about this software to pfc@stat.osu.edu.
##
## GNU GENERAL PUBLIC LICENSE, Version 3
## https://www.gnu.org/licenses/gpl-3.0.txt
## ======================================================================



## ===================================================================
## Purpose : Returns a vector of the times associated with the
##           scaling (is.wavelet=F) or wavelet (is.wavelet=T)
##           coefficients for level 'j' of the MODWT object 'dx'. 
## Updated : pfc@stat.osu.edu, Jun 2009.
## ===================================================================

modwt.times <- function (dx, j, is.wavelet=TRUE)
{
  start <- (sum(start(dx$x))-1)/frequency(dx$x)
  vj    <- dwt.filter.aligning(dx$filter, j, is.wavelet)
  rotate.vector(seq(from=start, by=deltat(dx$x), length=dx$N), abs(vj))
}




## ===================================================================
## Purpose : Returns the approximate bandpass frequencies for the
##           level 'js' wavelet coefficients of the modwt class 'dx'.
## Updated : pfc@stat.osu.edu, Apr 2005.
## ===================================================================

modwt.freqs <- function (dx, js)
{
  outer(c(0.5,1)/deltat(dx$x), 2^(-js))
}


