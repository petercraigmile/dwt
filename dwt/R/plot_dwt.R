## ======================================================================
## File     : dwt_plot.R
## Purpose  : R code for the Discrete Wavelet Transform (DWT).
## Requires : stackplot.R
## Updated  : pfc@stat.osu.edu, Mar 2009.
## 
## Copyright 2002--, Peter F. Craigmile, All Rights Reserved
## Address comments about this software to pfc@stat.osu.edu.
##
## GNU GENERAL PUBLIC LICENSE, Version 3
## https://www.gnu.org/licenses/gpl-3.0.txt
## ======================================================================


## ===================================================================
## Purpose : Produce a wavelet decomposition plot of the dwt
##           class 'dx' using a stacked plot.
## Updated : pfc@stat.osu.edu, Jun 2002.
## ===================================================================

plot.dwt <- function (x, xlab="time", col="black", bg="white",
                      bdcol="black", ...)
{
  xs <- ys <- vlines <- list()
  start      <- start(x$x)[1] + (start(x$x)[2]-1)/frequency(x$x)
  ts         <- seq(from=start, length=x$N, by=deltat(x$x))
  xlabels    <- NULL
  types      <- c("l", rep("h", x$nlevels), "l")
  xs[[1]]    <- ts
  ys[[1]]    <- x$x
  xlabels[1] <- "X"
  
  vn <- 0
  for (j in 1:x$nlevels)
  {
    ks           <- dwt.times(x, j, TRUE)
    xs[[j+1]]    <- sort(ks)
    ys[[j+1]]    <- x$W[[j]][order(ks)]
    xlabels[j+1] <- paste("W", j, sep="")

    if (x$Bj[j]>0)
    {
      delta <- diff(ks)[1]/2
      vlines[[vn+1]] <- c(j,ks[1]  - delta)
      vlines[[vn+2]] <- c(j,ks[x$Bj[j]] + delta)
      vn <- vn + 2
    }
  }

  ks                     <- dwt.times(x, x$nlevels, FALSE)
  xs[[x$nlevels+2]]      <- sort(ks)
  ys[[x$nlevels+2]]      <- x$V[order(ks)]
  xlabels[[x$nlevels+2]] <- paste("V",x$nlevels, sep="")

  stackplot(xs, ys, ts, xlabels, types, vlines, xlab=xlab,
            col=col, bg=bg, bdcol=bdcol, ...)
}

