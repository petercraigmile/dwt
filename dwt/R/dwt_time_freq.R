## ======================================================================
## File     : dwt_time_freq.R
## Purpose  : R code for aligning wavelet and scaling coefficients with time
##            or approximately with frequency in the DWT.
## Updated  : pfc@stat.osu.edu, Mar 2009.
## 
## Copyright 2002--, Peter F. Craigmile, All Rights Reserved
## Address comments about this software to pfc@stat.osu.edu.
##
## GNU GENERAL PUBLIC LICENSE, Version 3
## https://www.gnu.org/licenses/gpl-3.0.txt
## ======================================================================



## ===================================================================
## Purpose : Returns a vector of the times associated with the
##           scaling (is.wavelet=FALSE) or wavelet (is.wavelet=TRUE)
##           coefficients for level 'j' of the DWT object 'dx'. 
## Updated : pfc@stat.osu.edu, Jun 2002.
## ===================================================================

dwt.times <- function (dx, j, is.wavelet=TRUE)
{
  start <- start(dx$x)[1] + (start(dx$x)[2]-1)/frequency(dx$x)
  vj    <- dwt.filter.aligning(dx$filter, j, is.wavelet)
  start + deltat(dx$x) * ((2^j*(seq(dx$Nj[j])) - 1 - abs(vj)) %% dx$N)
}



## ===================================================================
## Purpose : Returns the approximate bandpass frequencies for the
##           level 'js' wavelet coefficients of the dwt class 'dx'.
## Updated : pfc@stat.osu.edu, Jun 2002.
## ===================================================================

dwt.freqs <- function (dx, js)
{
  outer(c(0.5,1)/deltat(dx$x), 2^(-js))
}




## ======================================================================
## Purpose : calculate the aligning factor for a level 'j' wavelet
##           (is.wavelet==TRUE) or scaling (is.wavelet==FALSE) coefficient
##           for a specific wavelet 'filter'.
## Updated : pfc@stat.osu.edu, Oct 2001. 
## ======================================================================

dwt.filter.aligning <- function (filter, j, is.wavelet=TRUE)
{
  L  <- filter$L
  Lj <- (2^j-1)*(L-1)+1
  L2 <- L/2
  ch1 <- substring(filter$name,1,1)
  
  ## Coiflets (Cx)
  if (ch1=="C") v <- -2*L/3+1
  ## Best Localised (BLx)
  else if (ch1=="B")
  {
    if (L==14)      v <- -5
    else if (L==18) v <- -11
    else if (L==20) v <- -9
  }
  ## Daubechies (Dx) and Least Asymmetric (LAx)
  else
  {
    if ((L == 4) | (L==6))  v <- -1
    else if ((L == 10) | (L==18)) v <- -L2
    else if (L == 14) v <- -L2 + 2
    else v <- -L2 + 1
  }

  ## equations 114(c) of Percival and Walden (2000)
  if (is.wavelet==TRUE) -(Lj/2 + L2 + v - 1)
  else               (Lj-1)*v/(L-1)
}


