## ======================================================================
## GNU GENERAL PUBLIC LICENSE, Version 3
## https://www.gnu.org/licenses/gpl-3.0.txt
## ======================================================================


plot.modwt <- function (x, xlab="time", col="black", 
		        bg="white", bdcol="blue", ...) {
  
## ===================================================================
## Purpose : to produce a wavelet decomposition plot of the MODWT
##           class 'dx' using a stacked plot.
## Updated : pfc@stat.osu.edu, April 2002.
## ===================================================================

  xs <- ys <- vlines <- list()
  ys         <- list()
  start      <- (sum(start(x$x))-1)/frequency(x$x)
  ts         <- seq(from=start, length=x$N, by=deltat(x$x))
  xlabels    <- rep("",x$nlevels+2)
  types      <- c("l", rep("l",x$nlevels), "l")
  xs[[1]]    <- ts
  ys[[1]]    <- x$x
  xlabels[1] <- "X"
  vn         <- 0

  for (j in 1:x$nlevels)
  {
    ks           <- modwt.times(x, j, TRUE)
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

  ks                     <- modwt.times(x, x$nlevels, FALSE)
  xs[[x$nlevels+2]]      <- sort(ks)
  ys[[x$nlevels+2]]      <- x$V[order(ks)]
  xlabels[[x$nlevels+2]] <- paste("V",x$nlevels, sep="")

  if (x$Bj[x$nlevels]>0)
  {
    delta <- diff(ks)[1]/2
    vlines[[vn+1]] <- c(x$nlevels+1,ks[1]  - delta)
    vlines[[vn+2]] <- c(x$nlevels+1,ks[x$Bj[x$nlevels]] + delta)
    vn <- vn + 2
  }
  
  stackplot(xs, ys, ts, xlabels, types, vlines, xlab=xlab,
            col=col, bg=bg, bdcol=bdcol, ...)
}

