
/* ======================================================================
// File     : dwt.c
// Purpose  : C code used to implement the DWT
//            along with their inverses
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


/* ======================================================================
// 'dwt.R' functions with circular boundary conditions
// ====================================================================*/

void R_dwt_cascade_previous (double *currV, double *currW, int *M, 
			     double *prevV, double *g, double *h, int *L)
{
  int s, u, t, l;
  double Vs, Vsp1;

  for (t = 0, s = 0; t < (*M); t++, s+=2)
  {
    u = t;    
    Vs = Vsp1 = 0.0;

    for (l = 0; l < (*L); l+=2)
    {
      Vs   += g[l+1]*currV[u] + h[l+1]*currW[u];
      Vsp1 += g[l]  *currV[u] + h[l]  *currW[u];

      u++;
      if (u>=(*M)) u = 0;
    }

    prevV[s]   = Vs;
    prevV[s+1] = Vsp1;
  }
}
 

void R_dwt_cascade_next (double *currV, int *M, 
			 double *nextV, double *nextW,
			 double *g, double *h, int *L)
{
  int t, u, l;
  double nextV_t, nextW_t;
  
  for (t=0; t<(*M)/2; t++)
  {
    u  = 2*t + 1;
    nextV_t = g[0] * currV[u];
    nextW_t = h[0] * currV[u];

    for (l=1;  l<(*L);  l++)
    {
      u--;
      if (u<0) u = (*M) - 1;
      nextV_t += g[l] * currV[u];
      nextW_t += h[l] * currV[u];
    }

    nextV[t] = nextV_t;
    nextW[t] = nextW_t;
  }
}





