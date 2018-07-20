//
// Created by rais on 23/09/17.
//

#ifndef GC_MATOPS_H
#define GC_MATOPS_H

#include <fftw3.h>
#include <algorithm>

void matConvolve2D(double *img, int w, int h, double *kernel, int wk, int hk, double *output, int borderPolicy=0);
void matMultiply(double *I1, int w1, int h1, double *I2, int w2, int h2, double **res);
inline double matDotProduct(double *a, double *b, int len, int step)
{
    double r = 0;
    while (len--) {
        r += *a++ * *b;
        b += step;
    }
    return r;
}
inline void complexMult(fftw_complex &n1, fftw_complex &n2, fftw_complex &res) {
    res[0] = n1[0] * n2[0] - n1[1] * n2[1];
    res[1] = n1[0] * n2[1] + n1[1] * n2[0];
}
inline void complexConj(fftw_complex &n, fftw_complex &res) {
    res[0] = n[0];
    res[1] = -n[1];
}

inline void matrix3by3Inverse(double a[9], double res[9]) {
    double determinant = 0;
    for(int i=0;i<3;i++)
        determinant += (a[i]*(a[3+((i+1)%3)]*a[6+((i+2)%3)] - a[3+((i+2)%3)]*a[6+((i+1)%3)]));
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++) {
            res[i + j * 3] = ((a[((i + 1) % 3) * 3 + ((j + 1) % 3)] * a[((i + 2) % 3) * 3 + ((j + 2) % 3)]) -
                              (a[((i + 1) % 3) * 3 + ((j + 2) % 3)] * a[((i + 2) % 3) * 3 + ((j + 1) % 3)])) /
                             determinant;
        }
    }
}
inline double complexMag(fftw_complex &n) {
    return sqrt(n[0] * n[0] + n[1] * n[1]);
}
inline void complexMag(fftw_complex *img, int w, int h, double *res, bool calcLog=true) {
    for (int x = 0; x < w; x++) {
        for (int y = 0; y < h; y++) {
            res[x + y * w] = complexMag(img[x+y*w]);
            if (calcLog) res[x + y * w] = log(res[x + y * w] + 1);
        }
    }
}


inline void circshift(double *out, const double *in, int w, int h, int xshift, int yshift)
{
    for (int i =0; i < w; i++) {
        int ii = (i + xshift) % w;
        for (int j = 0; j < h; j++) {
            int jj = (j + yshift) % h;
            out[ii * h + jj] = in[i * h + j];
        }
    }
}

inline void fftshift(double *out, const double *in, int w, int h) {
    circshift(out, in, w, h, (w/2), (h/2));
}

template<class T> void ifftShift(T *out, const T* in, size_t nx, size_t ny)
{
    const size_t hlen1 = (ny+1)/2;
    const size_t hlen2 = ny/2;
    const size_t shft1 = ((nx+1)/2)*ny + hlen1;
    const size_t shft2 = (nx/2)*ny + hlen2;

    const T* src = in;
    for(T* tgt = out; tgt < out + shft1 - hlen1; tgt += ny, src += ny) { // (nx+1)/2 times
        std::copy(src, src+hlen1, tgt + shft2);          //1->4
        std::copy(src+hlen1, src+ny, tgt+shft2-hlen2); } //2->3
    src = in;
    for(T* tgt = out; tgt < out + shft2 - hlen2; tgt += ny, src += ny ){ // nx/2 times
        std::copy(src+shft1, src+shft1+hlen2, tgt);         //4->1
        std::copy(src+shft1-hlen1, src+shft1, tgt+hlen2); } //3->2
};
#endif //GC_MATOPS_H
