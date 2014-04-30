## ======================================================================
## File     : rotate_vector.R
## Updated  : pfc@stat.osu.edu, Apr 2005.
## 
## Copyright 2002--, Peter F. Craigmile, All Rights Reserved
## Address comments about this software to pfc@stat.osu.edu.
##
## GNU GENERAL PUBLIC LICENSE, Version 3
## https://www.gnu.org/licenses/gpl-3.0.txt
## ======================================================================



## ======================================================================
## Purpose : Rotate the vector 'x' 'delta' units to the right.
## Updated : pfc@stat.osuedu, June 2002.
## ======================================================================

rotate.vector <- function (x, delta)
{
  n <- length(x)
  if (n > 0)
  {
    mod.delta <- n-(delta %% n)
    c(x,x)[(mod.delta+1):(mod.delta+n)]
  }
  else x
}

