//
// Created by rais on 23/09/17.
//

#ifndef GC_WINDOW_H
#define GC_WINDOW_H
void create1DTukeyWindow(double *res, int L, double cutoff);
void create2DTukeyWindow(double **res, int w, int h, double cutoff);
#endif //GC_WINDOW_H
