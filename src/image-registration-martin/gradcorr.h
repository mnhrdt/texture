//
// Created by rais on 23/09/17.
//

#ifndef GC_GRADCORR_H
#define GC_GRADCORR_H
void removeNans(double *I, int w, int h);
void height_shift(const char *I1, const char *I2, double resX, double resY, double *resZ);
bool registerGC(double *I1, double *I2, int w, int h, double *resX, double *resY);
bool registerGCFiles(const char *I1, const char *I2, double *resX, double *resY);

#endif //GC_GRADCORR_H
