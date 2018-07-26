//
// Created by rais on 19/09/17.
//

#ifndef GC_INTERPOLATION_H
#define GC_INTERPOLATION_H
double cubicInterpolate (double p[4], double x);

double bicubicInterpolate (double p[4][4], double x, double y);

double GetInterpolatedValueBicubic(double *img, int width, int height, double x, double y);

inline float BilinearInterpolation(float q11, float q12, float q21, float q22, float x1, float x2, float y1, float y2, float x, float y)
{
    double wx = (x-x1)/(x2-x1);
    double wy = (y-y1)/(y2-y1);
    double res =
            (1 - wy) * ((1-wx) * q11 + wx * q12) +
            (    wy) * ((1-wx) * q21 + wx * q22) ;

    return res;
}

inline float GetInterpolatedValueBilinear(double *img, int width, int height, double x, double y)
{
    (void)height;
    int xF = floor(x);
    int yF = floor(y);
    float res = BilinearInterpolation(
            img[xF + yF * width], img[xF+1 + yF * width],
            img[xF + (yF+1) * width], img[xF+1 + (yF+1) * width],
            xF, xF+1, yF, yF+1, x, y);
    return res;
}

#endif //GC_INTERPOLATION_H
