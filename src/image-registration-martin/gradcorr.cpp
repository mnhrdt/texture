//
// Created by rais on 23/09/17.
//
#include <cstdlib>
#include <fftw3.h>
#include <cstring>
#include <cmath>
#include <iostream>
#include <vector>
#include "gradcorr.h"
#include "gradient.h"
#include "window.h"
#include "matops.h"

extern "C" {
#include "iio.h"
}

void removeNans(double *I, int w, int h) {
    for (int x=0;x<w;x++) {
        for (int y=0;y<h;y++) {
            if (std::isnan(I[x+w*y]) || std::isinf(I[x+w*y]))
                I[x+w*y] = 0;
        }
    }
}

void crop_images(double *image1, double *image2, int w1, int w2, 
        double *img1, double *img2, int h)
{
    int w = std::min(w1,w2);

    for (int i = 0; i < w; i++)
        for (int j = 0; j < h; j++)
        {
            img1[i + j * w] = image1[i + j * w1];
            img2[i + j * w] = image2[i + j * w2];
        }
}

void median_height_estimation(double *image1, double *image2, double *imgs, 
        int w1, int w2, int w, int h, double resX, double resY, double *resZ){

    std::vector<double> tmp1;
    std::vector<double> tmp2;

    // loop over the points nonan in both image to fill the above lists
    // fill in registered image (height from image2)
    for (int i = 0; i < w; i++)
        for (int j = 0; j < h; j++){
            int i2 = i - resX;
            int j2 = j - resY;
            if (i2 < 0 || i2 > w - 1 || j2 < 0 || j2 > h - 1)
                continue;
            else 
                imgs[i+j*w] = image2[i2 + j2 * w2];
            if ( !std::isnan(image1[i+j*w1]) && !std::isnan(image2[i2+j2*w2]) ){
                tmp1.push_back(image1[i+j*w1]);
                tmp2.push_back(image2[i2+j2*w2]);
            }
        }

    // compute difference of medians
    std::sort(tmp1.begin(), tmp1.end());
    std::sort(tmp2.begin(), tmp2.end());
    resZ[0] = tmp2[tmp2.size()/2]-tmp1[tmp1.size()/2];

    // adjust height for registered image
    for (int i = 0; i < w*h; i++)
        imgs[i] -= resZ[0];

}

void mode_height_estimation(double *image1, double *imgs, 
        int w1, int w, int h, double *resZ,
        double *alpha){

    std::vector<int> hist(40,0);
    for (int i = 0; i < w; i++)
        for (int j = 0; j < h; j++)
            if (abs(imgs[i+j*w] - image1[i+j*w1]) < 4)
                for (int k = 0; k < 39; k++){
                    hist[k] += (imgs[i+j*w] - image1[i+j*w1] > - 3.9 + 0.2 * k) &&
                        (imgs[i+j*w] - image1[i+j*w1] < - 3.9 + 0.2 * (k + 1));
                }

    std::vector<int>::iterator result = std::max_element(hist.begin(), hist.end());
    int idx = std::distance(hist.begin(), result);

    alpha[0] = - 3.9 + 0.2 * (idx + 0.5);
    for (int i = 0; i < w*h; i++)
        imgs[i] -= alpha[0];

    resZ[0] += alpha[0];
    hist.clear();
}

void height_shift(const char *I1, const char *I2, double resX, double resY, double *resZ, 
        char *f_out){

    // read input images
    int w1, h1, pd1, w2, h2, pd2;
    double *image1 = (double *) iio_read_image_double_vec(I1, &w1, &h1, &pd1);
    double *image2 = (double *) iio_read_image_double_vec(I2, &w2, &h2, &pd2);
    
    // initialise registered image
    int w = std::min(w1,w2); int h = std::min(h1,h2);
    double *imgs = (double *) malloc(w * h * sizeof(double));

    median_height_estimation(image1, image2, imgs, w1, w2, w, h, resX, resY, resZ);
    double alpha = 10;

    mode_height_estimation(image1, imgs, w1, w, h, resZ, &alpha);

    // write registered image
    iio_write_image_double_vec(f_out, &(imgs[0]), w, h, 1);
    free(imgs);
}

bool registerGCFiles(const char *I1, const char *I2, double *resX, double *resY) {
    int w1, h1, pd1, w2, h2, pd2;
    double *image1 = iio_read_image_double_vec(I1, &w1, &h1, &pd1);
    double *image2 = iio_read_image_double_vec(I2, &w2, &h2, &pd2);
    
    int w = std::min(w1,w2); int h = std::min(h1,h2);
    double *img1 = (double *) malloc(w * h * sizeof(double));
    double *img2 = (double *) malloc(w * h * sizeof(double));
    crop_images(image1, image2, w1, w2, img1, img2, h); 
    
//    if (w1 != w2 || h1 != h2 || pd1 != pd2 || pd1 != 1) {
//        std::cout << "Images must have the same size, height and be single channel." << std::endl;
//        std::cout << "Image 1: (" << w1 << "x" << h1 << "x" << pd1 << std::endl;
//        std::cout << "Image 2: (" << w2 << "x" << h2 << "x" << pd2 << std::endl;
//    } else {
      return registerGC(img1, img2, w, h, resX, resY);
//    }
}

bool registerGC(double *I1, double *I2, int w, int h, double *resX, double *resY) {
    //double cutoffFreq = 0.25;
    // First, compute image gradients
    double *dx1 = (double *) malloc(w * h * sizeof(double));
    double *dy1 = (double *) malloc(w * h * sizeof(double));
    double *dx2 = (double *) malloc(w * h * sizeof(double));
    double *dy2 = (double *) malloc(w * h * sizeof(double));
    computeCenteredGradient(I1, w, h, dx1, dy1);
    computeCenteredGradient(I2, w, h, dx2, dy2);
    // Eliminate nans from gradients
    removeNans(dx1, w, h);
    removeNans(dx2, w, h);
    removeNans(dy1, w, h);
    removeNans(dy2, w, h);


    // Create complex gradient image
    fftw_complex *gI1, *gI2, *FFTgI1, *FFTgI2;
    gI1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * w * h);
    gI2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * w * h);
    FFTgI1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * w * h);
    FFTgI2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * w * h);
    for (int x=0; x<w; x++) {
        for (int y=0;y<h;y++) {
            gI1[y*w+x][0] = dx1[y*w+x];
            gI1[y*w+x][1] = dy1[y*w+x];
            gI2[y*w+x][0] = dx2[y*w+x];
            gI2[y*w+x][1] = dy2[y*w+x];
        }
    }
    // Create Tukey Window to avoid border effects on the FFT
    double *tukey;// = (double *) malloc(175 * sizeof(double));
    create2DTukeyWindow(&tukey, w, h, 0.25);
    // Apply tukey window in the spatial domain
    for (int x = 0; x < w; x++) {
        for (int y = 0; y < h; y++) {
            gI1[y*w+x][0] *= tukey[y*w+x];
            gI1[y*w+x][1] *= tukey[y*w+x];
            gI2[y*w+x][0] *= tukey[y*w+x];
            gI2[y*w+x][1] *= tukey[y*w+x];
        }
    }
    free(tukey);
    // Apply FFT
    fftw_plan p1, p2, p3;
    p1 = fftw_plan_dft_2d(h, w, gI1, FFTgI1, FFTW_FORWARD, FFTW_ESTIMATE);
    p2 = fftw_plan_dft_2d(h, w, gI2, FFTgI2, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p1);
    fftw_execute(p2);

    // Save magnitude
//    double *mag = (double *) malloc(w * h * sizeof(double));
//    char filename[200];
//    complexMag(FFTgI1, w, h, mag);
//    strcpy(filename, "/home/rais/PhD/Thales/data/maggI1.tiff");
//    iio_write_image_double(filename, mag, w, h);
//    complexMag(FFTgI2, w, h, mag);
//    strcpy(filename, "/home/rais/PhD/Thales/data/maggI2.tiff");
//    iio_write_image_double(filename, mag, w, h);
//    free(mag);

    // Now compute the correlation in the Fourier domain
    fftw_complex *FFTgc = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * w * h);
    fftw_complex *gcC = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * w * h);
    fftw_complex conjgI2;
    for (int x = 0; x < w; x++) {
        for (int y = 0; y < h; y++) {
            complexConj(FFTgI2[x + w * y], conjgI2);
            complexMult(FFTgI1[x + w * y], conjgI2, FFTgc[x + w * y]);
        }
    }
    p3 = fftw_plan_dft_2d(h, w, FFTgc, gcC, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p3);
    double *gc = (double *) malloc(w * h * sizeof(double));
    for (int x = 0; x < w; x++) {
        for (int y = 0; y < h; y++) {
            gc[x + w * y] = gcC[x + w * y][0];
        }
    }
    //strcpy(filename, "/home/rais/PhD/Thales/data/GC.tiff");
    //iio_write_image_double(filename, gc, w, h);

    double *gcShifted = (double *) malloc(w * h * sizeof(double));
    ifftShift<double>(gcShifted, gc, h, w);
//    fftshift(gcShifted, gc, w, h);
    free(gc);
    gc = gcShifted;
    //fftshift(gc, w, h);
    //strcpy(filename, "/home/rais/PhD/Thales/data/GC2.tiff");
    //iio_write_image_double(filename, gc, w, h);

    // Find the peak
    double max = 0;
    int mx, my;
    for (int x = 0; x < w; x++) {
        for (int y = 0; y < h; y++) {
            if (fabs(gc[x + w*y]) > max) {
                max = fabs(gc[x + w * y]);
                mx = x;
                my = y;
            }
        }
    }
    // Interpolate the peak to achieve better accuracy
    bool failed = false;
    // Avoid peaks on the image borders
    if (mx < w - 2 && my < h-2 && mx > 0 && my > 0) {
        double sy = (0.5) * (gc[mx + (my + 1) * w] - gc[mx + (my - 1) * w]) /
                    (2.0 * gc[mx + my * w] - gc[mx + (my - 1) * w] - gc[mx + (my + 1) * w]);
        double sx = (0.5) * (gc[mx + 1 + my * w] - gc[mx - 1 + my * w]) /
                    (2.0 * gc[mx + my * w] - gc[mx + 1 + my * w] - gc[mx - 1 + my * w]);
        double resy = (double) my + sy;
        double resx = (double) mx + sx;
        *resX = floor((double) w/2.0) - resx;
        *resY = floor((double) h/2.0) - resy;
    }
    else {
        failed = true;
    }

    fftw_destroy_plan(p1);
    fftw_destroy_plan(p2);
    fftw_destroy_plan(p3);
    fftw_free(FFTgc);
    fftw_free(gcC);
    fftw_free(gI1);
    fftw_free(gI2);
    fftw_free(FFTgI1);
    fftw_free(FFTgI2);
    return !failed;

}
