## ======================================================================
## File    : dwt_filters_cascades.R
## Purpose : R code to define and manipulate DWT scaling filters
##           at many different wavelet levels.
## Updated : pfc@stat.osu.edu
## 
## Copyright 2002--2009, Peter F. Craigmile, All Rights Reserved
## Address comments about this software to pfc@stat.osu.edu.
##
## GNU GENERAL PUBLIC LICENSE, Version 3
## https://www.gnu.org/licenses/gpl-3.0.txt
## ======================================================================



## ======================================================================
## Purpose : Using the dwt.filter object 'filter', calculate the
##           all the wavelet and scaling filters up to level 'J'.
## Assumes : 'J' is an integer > 0.
## Updated : pfc@stat.osu.edu, Mar 2003.
## ======================================================================

dwt.all.filters <- function (filter, J) {

  filter.list <- list(filter)
  
  if (J > 1) 
    for (j in 2:J) {
      
      filter.list[[j]] <- dwt.next.filters(filter.list[[j-1]], filter, j)
    }
  
  filter.list
}




## ======================================================================
## Purpose : Given the wavelet and scaling filters at level 'j-1'
##           'prev.filters' calculate the level 'j' filters
##           using the dwt.filter object 'filter'.
## Assumes : 'j' is an integer > 1.
## Updated : pfc@stat.osu.edu, Jun 2010.
## =====================================================================

dwt.next.filters <- function (prev.filters, filter, j) {
  
  Lj <- (2^j-1) * (filter$L - 1) + 1
  
  next.scaling <- next.wavelet <- rep(NA, Lj)
  
  ks <- seq(0, length(prev.filters$scaling)-1)
  
  for (l in 1:Lj) {
    
    ms <- l-2*ks
    sel <- ms>=1 & ms<=filter$L
    
    next.scaling[l] <- sum(filter$scaling[ms[sel]] * prev.filters$scaling[sel])
    next.wavelet[l] <- sum(filter$scaling[ms[sel]] * prev.filters$wavelet[sel])
  }
  
  list(scaling=next.scaling,
       wavelet=next.wavelet)
}
  
