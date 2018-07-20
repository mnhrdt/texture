//
// Created by rais on 23/09/17.
//

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include "window.h"
#include "matops.h"

void create1DTukeyWindow(double *res, int n, double cutoff)
{
   // Defines period of the taper as 1/2 period of a sine wave.
    double per = cutoff/2;
    int tl = floor(per*(n-1))+1;
    int th = n-tl+1;
    // Window is defined in three sections: taper, constant, taper
    for (int i=0;i<n;i++) {
        double pos = (double) i / (n-1);
        if (i < tl) {
            res[i] = ((1+cos(M_PI /per * (pos - per)))/2);
        } else if (i < th - 1) {
            res[i] = 1;
        } else {
            res[i] = ((1+cos(M_PI / per * (pos - 1 + per)))/2);
        }
    }
}

void create2DTukeyWindow(double **res, int w, int h, double cutoff)
{
    double *wc = (double *) malloc(h * sizeof(double));
    double *wr = (double *) malloc(w * sizeof(double));
    create1DTukeyWindow(wc, h, cutoff);
    create1DTukeyWindow(wr, w, cutoff);
    matMultiply(wc, 1, h, wr, w, 1, res);
    free(wc);
    free(wr);
}