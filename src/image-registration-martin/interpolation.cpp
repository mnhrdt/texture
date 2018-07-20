#include <cmath>
#include "interpolation.h"

double cubicInterpolate (double p[4], double x) {
    return p[1] + 0.5 * x*(p[2] - p[0] + x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + x*(3.0*(p[1] - p[2]) + p[3] - p[0])));
}

double bicubicInterpolate (double p[4][4], double x, double y) {
    double arr[4];
    arr[0] = cubicInterpolate(p[0], y);
    arr[1] = cubicInterpolate(p[1], y);
    arr[2] = cubicInterpolate(p[2], y);
    arr[3] = cubicInterpolate(p[3], y);
    return cubicInterpolate(arr, x);
}

double GetInterpolatedValueBicubic(double *img, int width, int height, double x, double y) {
    double arr[4][4];
    int i = floor(x);
    int j = floor(y);
    for (int l = -2; l < 2; l++) {
        for (int k = -2; k < 2; k++) {
            int pos = (i + l) + (j + k) * width;
            if (pos >= 0 && pos <= width*height)
                arr[l+2][k+2] = img[(i + l) + (j + k) * width];
        }
    }
    double res = bicubicInterpolate(arr, x / (double) width, y / (double) height);
    return res;
}


