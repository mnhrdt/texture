//
// Created by rais on 29/09/17.
//

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <cmath>
#include "satreg.h"
extern "C" {
#include "iio.h"
}
#include "gradcorr.h"
#include "gradient.h"
#include "interpolation.h"
#include "matops.h"


void register_satellite_aerial(const char *satImg, const char *aerImg, const char *opt2ref, double &deltaX, double &deltaY,
                               double **res)
{
    int w, h, pd;
    double *img = iio_read_image_double_vec(satImg, &w, &h, &pd);
    std::cout << "Satellite image  size: " << w << "x" << h << "x" << pd << std::endl;
    double *imgGray = (double*) malloc(w * h * sizeof(double));
    // Compute the mean of the channels
    for (int i=0;i<w;i++) {
        for (int j=0;j<h;j++) {
            double acum = 0.0;
            for (int k=0;k<pd;k++) {
                acum += img[(i+j*w)*pd+k];
            }
            imgGray[i + j*w] = acum / 3.0;
        }
    }
    int wOpt, hOpt;
    double *imgOpt = iio_read_image_double(aerImg, &wOpt, &hOpt);
    std::cout << "Aerial image size: " << wOpt << "x" << hOpt << std::endl;

    // Parse homography from file
    double H[9];
    parseHomography(H, opt2ref);
    double Hinv[9];
    matrix3by3Inverse(H, Hinv);
    std::cout << "H inv: ";
    for (int i=0;i<9;i++) {
        Hinv[i] /= Hinv[8];
        std::cout << Hinv[i] << " ";
    }
    std::cout << std::endl;


    int x1, x2, y1, y2;
    double *imOptResampled = ResampleImageHomo(imgOpt, wOpt, hOpt, w, h, Hinv, &x1, &x2, &y1, &y2);
//    iio_write_image_double("imOptResampled.tif", imOptResampled, w, h);

    // Register images using only rectangle (x1,y1, x2, y2)
    int szX = x2 - x1 + 1;
    int szY = y2 - y1 + 1;
    double *imBase = (double *) malloc(szX * szY * sizeof(double));
    double *imToSearch = (double *) malloc(szX * szY * sizeof(double));
    for (int x = 0; x < szX; x++) {
        for (int y =0; y < szY; y++) {
            imBase[x + y*szX] = imgGray[x+x1+(y+y1)*w];
            imToSearch[x + y*szX] = imOptResampled[x+x1+(y+y1)*w];
        }
    }
//    iio_write_image_double("imBase.tif", imBase, szX, szY);
//    iio_write_image_double("imToSeach.tif", imToSearch, szX, szY);

    // Compute gradients from both images
    double *gradMagImToSearch = (double *) malloc(szX * szY * sizeof(double));
    double *gradMagImBase = (double *) malloc(szX * szY * sizeof(double));
    sobelMagnitude(imToSearch, szX, szY, gradMagImToSearch);
    sobelMagnitude(imBase, szX, szY, gradMagImBase);
    for (int x=0;x<szX;x++) {
        for (int y=0;y<szY;y++) {
            if (std::isnan(gradMagImToSearch[x+y*szX]))
                gradMagImToSearch[x+y*szX] = 0;
        }
    }

//    iio_write_image_double("gradMagImToSearch.tif", gradMagImToSearch, szX, szY);
//    iio_write_image_double("gradMagImBase.tif", gradMagImBase, szX, szY);

    registerGC(gradMagImToSearch, gradMagImBase, szX, szY, &deltaX, &deltaY);

    std::cout << "Estimated displacement: (" << deltaX << "," << deltaY << ")" << std::endl << std::endl;

    if (res != NULL) {
        H[2] += deltaX;
        H[5] += deltaY;
        matrix3by3Inverse(H, Hinv);
        *res = ResampleImageHomo(imgOpt, wOpt, hOpt, w, h, Hinv, &x1, &x2, &y1, &y2);
//        double H2[9] = {1, 0, -deltaX, 0, 1, -deltaY, 0, 0, 1};
//        *res = ResampleImageHomo(imOptResampled, w, h, w, h, H2, &x1, &x2, &y1, &y2);
    }

    free(gradMagImToSearch);
    free(gradMagImBase);
    free(imOptResampled);
    free(imBase);
    free(imToSearch);
    free(imgOpt);
    free(img);
    free(imgGray);

}

double *ResampleImageHomo(double *img, int wImg, int hImg, int wDst, int hDst, double H[9], int *minX,
                          int *maxX, int *minY, int *maxY) {
    double *res = (double *) malloc(wDst * hDst * sizeof(double));
    *minX = wDst;
    *minY = hDst;
    *maxX = *maxY = 0;
    int d = 1;
    for (int x = d; x < wDst + d; x++) {
        for (int y = d; y < hDst + d; y++) {
            double projZ = H[6] * (double) x + H[7] * (double) y + H[8];
            double projX = (H[0] * (double) x + H[1] * (double) y + H[2]) / projZ;
            double projY = (H[3] * (double) x + H[4] * (double) y + H[5]) / projZ;
            if (projX >= d && projX <= wImg + d - 1 && projY >= d && projY <= hImg + d - 1) {
                if (x < *minX) *minX = x;
                if (x > *maxX) *maxX = x;
                if (y < *minY) *minY = y;
                if (y > *maxY) *maxY = y;
                res[(x - d) + (y - d) * wDst] = GetInterpolatedValueBilinear(img, wImg, hImg, projX - d, projY - d);
                //res[(x-d) + (y-d) * wDst] = GetInterpolatedValueBicubic(img, wImg, hImg, projX, projY);
            } else {
                res[(x - d) + (y - d) * wDst] = nan("");
            }
        }
    }
    if (d != 0) {
        *minX -= d;
        *minY -= d;
        *maxX -= d;
        *maxY -= d;
    }
    return res;
}

const char* getfield(char* line, int num)
{
    const char* tok;
    for (tok = strtok(line, ","); tok && *tok; tok = strtok(NULL, ",\n")) {
        if (!--num)
            return tok;
    }
    return NULL;
}

void parseHomography(double H[9], const char *opt2refFile) {

    FILE *f = fopen(opt2refFile, "rt");
    char line[1024];
    char *token;
    int i = 0;
    while (fgets(line, 1024, f))
    {
        if (line[0] !=  '\n') {
            for (int k = 0; k < 3; k++) {
                char *tmp = strdup(line);
                const char *number = getfield(tmp, k + 1);
                double n = atof(number);
                H[3 * i + k] = n;
                free(tmp);
            }
        }
        i++;
    }
    fclose(f);
}

