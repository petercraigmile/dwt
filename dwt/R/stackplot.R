## ======================================================================
## Purpose: stackplot.R
## History: pfc@stat.osu.edu, Jun 2002.
##
## GNU GENERAL PUBLIC LICENSE, Version 3
## https://www.gnu.org/licenses/gpl-3.0.txt
## ======================================================================

## ======================================================================

## ======================================================================
## Purpose: [THIS DOCUMENTATION NEEDS TO BE FILLED IN]
## History: pfc@stat.osu.edu, May 2002.
## ======================================================================

stackplot <- function (xs, ys, ts, xlabels,
                       types, vlines=list(),
                       lowers=sapply(ys,min),
                       uppers=sapply(ys,max),
                       xlab="time",
                       col="black", bg="white", bdcol="blue", ...)
{
  par(bg=bg, col=col)

  mids    <- 0.5*(lowers+uppers)
  ranges  <- uppers-lowers
  cranges <- c(0,cumsum(ranges))
  n       <- length(ranges)

  plot(ts, seq(0,cranges[n+1],length=length(ts)),
       xlab=xlab, ylab="", type="n", yaxt="n",
       col=col, col.axis=col, col.lab=col, bg=bg)

  for (k in 1:n)
  {
    by <- cranges[k]-lowers[k]
    mtext(xlabels[k], side=2, at=by+mids[k], line=1, srt=90, col=col, cex=0.7)

    if (types[k]=="l")
      lines(xs[[k]],by+ys[[k]],col=col)
    else if (types[k]=="p")
      points(xs[[k]],by+ys[[k]],col=col)
    else if (types[k]=="h")
      segments(xs[[k]],by+ys[[k]],xs[[k]],by,col=col)

    abline(h=cranges[k],lty=1,col="grey")
  }
  abline(h=cranges[n+1],lty=1,col="grey")

  if (length(vlines)>0)
    for (k in 1:length(vlines))
    {
      j  <- vlines[[k]][1]
      x  <- vlines[[k]][2]
      by <- cranges[j+1]
      ty <- by + ranges[j+1]
      segments(x, by, x, ty, col=bdcol, ...)
    }
  
  invisible()
}


