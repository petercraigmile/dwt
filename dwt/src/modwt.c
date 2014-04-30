
/* ======================================================================
// File     : modwt.c
// Purpose  : C code used to implement MODWT
//            aint with their inverses
// Note     : Based on pseudo code from Percival and Walden, 
//            "Wavelet Methods for Time Series Analysis", (2000).
// Updated  : pfc@stat.osu.edu, Jun 2009
// 
// Copyright 2002--2009, Peter F. Craigmile, All Rights Reserved
// Address comments about this software to pfc@stat.osu.edu.
//
// GNU GENERAL PUBLIC LICENSE, Version 3
// https://www.gnu.org/licenses/gpl-3.0.txt
// ======================================================================*/

#include <R.h>
#include <Rmath.h>



void R_modwt_cascade_next (double *currV, int *N, int *j,
			   double *nextV, double *nextW,
			   double *g, double *h, int *L)
{
  int t, k, l, tau;
  double nextV_t, nextW_t;

  tau = 1 << (*j-1);

  for (t=0; t<(*N); t++)
  {
    k = t;
    nextV_t = g[0] * currV[k];
    nextW_t = h[0] * currV[k];

    for (l=1; l<(*L); l++)
    {
      k -= tau;
      while (k<0) k+=(*N);
      nextV_t += g[l] * currV[k];
      nextW_t += h[l] * currV[k];
    }

    nextV[t] = nextV_t;
    nextW[t] = nextW_t;
  }
}



void R_modwt_cascade_previous (double *currV, double *currW, 
			       int *N, int *tau, double *prevV,
			       double *g, double *h, int *L)
{
  int t, k, l;
  double prevV_t;

  for (t=0; t<(*N); t++)
  {
    k = t;
    prevV_t = g[0] * currV[k] + h[0] * currW[k];

    for (l=1; l<(*L); l++)
    {
      k += (*tau);
      while (k>=(*N)) k -= (*N);
      prevV_t += g[l] * currV[k] + h[l] * currW[k];
    }

    prevV[t] = prevV_t;
  }
}






