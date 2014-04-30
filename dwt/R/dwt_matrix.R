## ======================================================================
## Copyright 1998--2010, Peter F. Craigmile, All Rights Reserved
## Address comments about this software to pfc@stat.osu.edu.
##
## GNU GENERAL PUBLIC LICENSE, Version 3
## https://www.gnu.org/licenses/gpl-3.0.txt
## ======================================================================


## ======================================================================
## File     : dwt_matrix.R
## Contains : Calculate the orthogonal wavelet matrices
## Version  : 1.01
## Updated  : pfc@stat.osu.edu, Apr 2010.
## Requires : 'dwt.R', 'dwt_modwt_cascades.R', 'dwt_filters.R' and 'filters.R'
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



dwt.matrix <- function (dx, use.C=TRUE) {
  ## ======================================================================
  ## Calculates the dx$N x dx$N orthonormal dwt matrix using dwt class 'dx'.
  ## ADVICE: This matrix is sparse -- do not use for computer intensive
  ##         calculations.
  ## ======================================================================

  dz <- dwt(rep(0,dx$N), dx$filter$name, dx$nlevels, use.C=use.C)
  
  W <- matrix(0, dz$N, dz$N)
  r <- 1
  for (j in 1:dz$nlevels)
    for (k in 1:dz$Nj[j]) {
      dz$W[[j]][k] <- 1
      W[r, ] <- idwt(dz, use.C=use.C)
      dz$W[[j]][k] <- 0
      r <- r + 1
    }
  for (k in 1:dz$Nj[dz$nlevels]) {
    dz$V[k] <- 1
    W[r, ] <- idwt(dz, use.C=use.C)
    dz$V[k] <- 0
    r <- r + 1
  }
  
  W
}



dwt.matrix.row.subset <- function (dx, rows, use.C=TRUE) {
  ## ======================================================================
  ## Calculates the length(rows) x dx$N orthonormal dwt matrix using
  ## dwt class 'dx' at a subset of rows as given by 'rows'.
  ## ADVICE: This matrix is sparse -- do not use for computer intensive
  ##         calculations.
  ## ======================================================================

  dz <- dwt(rep(0,dx$N), dx$filter$name, dx$nlevels, use.C=use.C)
  
  sorted.rows <- sort(rows)
  W <- matrix(0, length(rows), dz$N)
  
  r <- 1
  s <- 1
  for (j in 1:dz$nlevels) for (k in 1:dz$Nj[j]) {
    
    if (r==sorted.rows[s]) {
      dz$W[[j]][k] <- 1
      W[s, ] <- idwt(dz, use.C=use.C)
      dz$W[[j]][k] <- 0
      s <- s + 1
    }      
    r <- r + 1
  }
  for (k in 1:dz$Nj[dz$nlevels]) {
    
    if (r==sorted.rows[s]) {
      dz$V[k] <- 1
      W[s, ] <- idwt(dz, use.C=use.C)
      dz$V[k] <- 0
      s <- s + 1
    }
    r <- r + 1
  }
  W
}




dwt.matrix.elems <- function (dx, use.C=TRUE) {
  ## ======================================================================
  ## Purpose : Calculate a subset of wavelet transform matrix for the dwt
  ##           object 'dx'.  Only include the rows corresponding
  ##           to the first wavelet coefficient per level, and the first
  ##           level dx$nlevels scaling coefficient.
  ## Updated : pfc@stat.osu.edu, Oct 2001. 
  ## ======================================================================

  dz <- dwt(rep(0,dx$N), dx$filter$name, dx$nlevels, use.C=use.C)
  W  <- matrix(0, dx$nlevels+1, dx$N)

  for (j in 1:dx$nlevels)
  {
    dz$W[[j]][1] <- 1
    W[j,]        <- idwt(dz)
    dz$W[[j]][1] <- 0
  }

  dz$V[1]    <- 1
  W[dx$nlevels+1,] <- idwt(dz)
  W	
}


