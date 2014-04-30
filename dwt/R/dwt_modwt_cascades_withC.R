## ======================================================================
## Copyright 1998--2010, Peter F. Craigmile, All Rights Reserved
## Address comments about this software to pfc@stat.osu.edu.
##
## GNU GENERAL PUBLIC LICENSE, Version 3
## https://www.gnu.org/licenses/gpl-3.0.txt
## ======================================================================


## ======================================================================
## File     : dwt_modwt_cascades.R
## Contains : Defines the filter cascades used in the DWT and MODWT.
## Version  : 1.01
## Updated  : pfc@stat.osu.edu, Apr 2010.
## Requires : 'dwt_filters.R' and 'filters.R'
##
## Notes    : This R code does not require extra dynamically allocated
##            C code, and so is very portable, but several times slower
##            than dedicated C code would be.
##            I could have optimized the code more, but I left the
##            filtering details intact for educational purposes.
##
## Reference:
##
##     D. B. Percival and A. T. Walden (2000),
##     Wavelet Methods for Time Series Analysis.
##     Cambridge, England: Cambridge University Press.
## ======================================================================



dwt.cascade.next <- function (V, dwt.filter.obj, use.C=TRUE) {
  ## ======================================================================
  ## Given the scaling coefficients 'V' on level j-1,
  ## calculate the scaling '$V' and wavelet coefficients '$W'
  ## on the level j using the dwt filter object 'dwt.filter.obj'.
  ## ======================================================================

  n <- length(V)
  if (n%%2!=0)
    stop("The vector V must be of even length")

   if (use.C) {

    if (!is.loaded("R_dwt_cascade_next", PACKAGE="dwt"))
      stop("The C function 'R_dwt_cascade_next' must be dynamically loaded.")
    
    z <- double(floor(n/2))
    results <- .C("R_dwt_cascade_next",
                  as.double(V), as.integer(n),
                  V = as.double(z),
                  W = as.double(z),
                  as.double(dwt.filter.obj$scaling),
                  as.double(dwt.filter.obj$wavelet),
                  as.integer(dwt.filter.obj$L),
                  PACKAGE="dwt")

    list(V=results$V, W=results$W)    
  } else {
    
    ## Optimization step: 'sel' is the vector that does the downsampling
    ## (I could also have used the function 'downsample.by.two').
    sel <- seq(2, n, 2)
    next.V <- circularly.filter(V, dwt.filter.obj$scaling, dwt.filter.obj$L)[sel]
    next.W <- circularly.filter(V, dwt.filter.obj$wavelet, dwt.filter.obj$L)[sel]
    
    list(V=next.V, W=next.W)
  }  
}




dwt.cascade.previous <- function (V, W, dwt.filter.obj, use.C=TRUE) {
  ## ======================================================================
  ## Given the scaling coefficients 'V' and wavelet coefficients 'W'
  ## on scale j+1, calculate the level j scaling coefficients using
  ## the dwt filter object 'dwt.filter.obj'.
  ## ======================================================================

  if (use.C) {

    if (!is.loaded("R_dwt_cascade_previous", PACKAGE="dwt"))
    stop("The C function 'R_dwt_cascade_previous' must be dynamically loaded.")
    
    Nj <- length(W)
    z <- double(Nj * 2)
    .C("R_dwt_cascade_previous",
       as.double(V), as.double(W), as.integer(Nj), V = z,
       as.double(dwt.filter.obj$scaling),
       as.double(dwt.filter.obj$wavelet),
       as.integer(dwt.filter.obj$L),
       PACKAGE="dwt")$V
    
  } else {
    V.up <- rev(upsample.by.two(V))
    W.up <- rev(upsample.by.two(W))
    
    rev(circularly.filter(V.up, dwt.filter.obj$scaling) +
        circularly.filter(W.up, dwt.filter.obj$wavelet))
  }
}


modwt.cascade.next <- function (V, j, dwt.filter.obj, use.C=TRUE) {
  ## ======================================================================
  ## Given the MODWT scaling coefficients 'V' on level 'j'-1,
  ## calculate the MODWT scaling '$V' and wavelet coefficients '$W'
  ## on the level 'j' using the dwt filter object 'dwt.filter.obj'.
  ## NOTE: Not been optimized for speed yet!
  ## ======================================================================

  N <- length(V)
  g <- dwt.filter.obj$scaling/sqrt(2)
  h <- dwt.filter.obj$wavelet/sqrt(2)

  if (use.C) {
    z <- double(N)

    results <- .C("R_modwt_cascade_next",
                  as.double(V), as.integer(N), as.integer(j), V = z, W = z,
                  as.double(g),
                  as.double(h),
                  as.integer(dwt.filter.obj$L),
                  PACKAGE="dwt")

    list(V=results$V, W=results$W)
  }
  else {
    tau <- 2^(j-1) - 1
    
    next.V <- circularly.filter(V,
                                periodize.filter(upsample.by(g, tau), N))
    next.W <- circularly.filter(V,
                                periodize.filter(upsample.by(h, tau), N))
    
    list(V=next.V, W=next.W)
  }
}


modwt.cascade.previous <- function (V, W, j, dwt.filter.obj, use.C) {
  ## ======================================================================
  ## Given the scaling coefficients 'V' and wavelet coefficients 'W'
  ## on scale j+1, calculate the level j scaling coefficients using
  ## the modwt filter object 'dwt.filter.obj'.
  ## PROBLEM: Really slow -- needs to be optimized!
  ## ======================================================================

  N <- length(V)
  g <- dwt.filter.obj$scaling/sqrt(2)
  h <- dwt.filter.obj$wavelet/sqrt(2)
  tau <- 2^(j-1)

  if (use.C) {
    z <- double(N)

    .C("R_modwt_cascade_previous",
       as.double(V), as.double(W), as.integer(N), as.integer(tau),
       prev.V=z,
       as.double(g),
       as.double(h),
       as.integer(dwt.filter.obj$L),
       PACKAGE="dwt")$prev.V
  }
  else {
    ls <- tau * (0 : (length(g)-1))
    
    sapply(0:(N-1), function (t) {
      ks <- ((t + ls) %% N) + 1
      sum(g * V[ks] + h * W[ks])
    })
  }
}


