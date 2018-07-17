//
// Created by rais on 23/09/17.
//

#include <cstdlib>
#include <cmath>
#include <cstring>
#include "matops.h"

extern "C" {
#include "iio.h"
}

void sobelMagnitude(double *img, int w, int h, double *output) {
    // Compute vertical and horizontal gradients
    double *Gx = (double *) malloc(w * h * sizeof(double));
    double *Gy = (double *) malloc(w * h * sizeof(double));
    double kx[9] = {1.0, 0.0, -1.0, 2.0, 0.0, -2.0, 1.0, 0.0, -1.0};
    double ky[9] = {1.0, 2.0, 1.0, 0.0, 0, 0.0, -1.0, -2.0, -1.0};
    matConvolve2D(img, w, h, kx, 3, 3, Gx, 0);
    matConvolve2D(img, w, h, ky, 3, 3, Gy, 0);
    for (int i=0;i<w*h;i++)
        output[i] = sqrt(Gx[i]*Gx[i]+Gy[i]*Gy[i]);
    free(Gx);
    free(Gy);
}

void computeCenteredGradient(double *img, int w, int h, double *dx, double *dy) {
    double k[3] = {1.0, 0.0, -1.0};
    matConvolve2D(img, w, h, k, 1, 3, dx, 0);
    matConvolve2D(img, w, h, k, 3, 1, dy, 0);
}
