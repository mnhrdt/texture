//
// Created by rais on 23/09/17.
//
#ifndef GC_GRADIENT_H
#define GC_GRADIENT_H

void sobelMagnitude(double *img, int W, int H, double *output);

void computeCenteredGradient(double *img, int W, int H, double *dx, double *dy);

#endif //GC_GRADIENT_H
