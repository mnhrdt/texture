//
// Created by rais on 29/09/17.
//

#ifndef GC_SATREG_H
#define GC_SATREG_H

void register_satellite_aerial(const char *satImg, const char *aerImg, const char *opt2ref, double &deltaX, double &deltaY,
                               double **res=NULL);

double *ResampleImageHomo(double *img, int wImg, int hImg, int wDst, int hDst, double H[9], int *minX,
                          int *maxX, int *minY, int *maxY);

void parseHomography(double H[9], const char *opt2refFile);
#endif //GC_SATREG_H
