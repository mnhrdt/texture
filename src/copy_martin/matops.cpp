//
// Created by rais on 23/09/17.
//

#include <stdlib.h>
#include <cmath>
#include "matops.h"

void matConvolve2D(double *img, int w, int h, double *kernel, int wk, int hk, double *output, int borderPolicy) {
    int posX, posY;
    for (int row = 0; row < h; row++) {
        for (int col = 0; col < w; col++) {
            double accumulation = 0;
            //double weightsum = 0;
            bool invalid = false;
            int i0 = (int) (wk / 2.0);
            int j0 = (int) (hk / 2.0);
            for (int i = -i0; i <= i0 && !invalid; i++) {
                for (int j = -j0; j <= j0 && !invalid; j++) {
                    posX = col + j;
                    posY = row + i;
                    switch (borderPolicy) {
                        // Replicate
                        case 0:
                            if (posX < 0) posX = 0;
                            if (posX >= w) posX = w-1;
                            if (posY < 0) posY = 0;
                            if (posY >= h) posY = h-1;
                            break;
                        // Symmetric
                        case 1:
                            if (posX < 0) posX = abs(posX);
                            if (posX >= w) posX -= (posX - w + 2);
                            if (posY < 0) posY = abs(posY);
                            if (posY >= h) posY -= (posY - h + 2);
                            break;
                    }

                    double k = img[posY * w + posX];
                    if (std::isnan(k)) {
                        invalid = true;
                        break;
                    }

                    // Check for 1D kernel
                    int realWK = wk;
                    if (hk == 1)
                        realWK = 1;
                    double kv = kernel[(i0 + i) * realWK + j0 + j];
                    accumulation += k * kv;
                    //weightsum += kv;
                }
            }
            if (!invalid)
                output[row * w + col] = accumulation; // / weightsum;
            else
                output[row * w + col] = nan(""); // / weightsum;
        }
    }
}


void matMultiply(double *I1, int w1, int h1, double *I2, int w2, int h2, double **res) {
    double *p, *pa;
    int i, j;
    if (w1 != h2) *res = NULL;
    *res = (double *) malloc(w2 * h1 * sizeof(double));
    p = *res;
    for (pa = I1, i = 0; i < h1; i++, pa += w1)
        for (j = 0; j < w2; j++) {
            *p++ = matDotProduct(pa, I2 + j, w1, w2);
        }
}
