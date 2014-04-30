
dwt.as.vector <- function (dx) {
  ## ======================================================================
  ## Take the wavelet and scaling coefficients in dwt object 'dx',
  ## and put them in a vector ordered as:
  ##
  ##   (W_1, W_2, ..., W_nlevels, V_nlevels)^T
  ## ======================================================================

  c(unlist(dx$W), dx$V)
}


dwt.map.vector <- function (dx, w) {
  ## ======================================================================
  ## Takes a vector of the form
  ##
  ##    (W_1, W_2, ..., W_nlevels, V_nlevels)^T
  ##
  ## and places these wavelet and scaling coefficients back in the
  ## dwt object 'dx'.
  ## ======================================================================

  z <- w
  
  for (j in 1:dx$nlevels) {

    sel <- 1:dx$Nj[j]
    dx$W[[j]] <- z[sel]
    z <- z[-sel]
  }

  dx$V <- z
  dx
}
