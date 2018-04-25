#ifndef COLORIZE_VERTICES_FROM_ONE_IMAGE_H
#define COLORIZE_VERTICES_FROM_ONE_IMAGE_H

#include <gdal.h>

#include "rpc.h"
#include "trimesh.h"

// ??????????????????????? pourquoi int pour z 
void lonlat_from_eastnorthzone(double out_lonlat[2], double e, double n, int z);

// cubic interpolation in dimension 1
static float cubic_interpolation(float v[4], float x);

void colorize(
        double *out_pan,               // vertices intensity
        double *out_msi,               // vertices multispectral infos
        double *out_rgb,               // vertices rgb colors
        struct trimesh *m,             // scaled mesh
        double *coord,                 // vertices utm coordinates
        struct rpc *pan_rpc,           // pan rpc
        struct rpc *msi_rpc,           // msi rpc
        int signed_zone,               // utm signed zone
        GDALRasterBandH *huge_pan_img, // satellite pan image 
        GDALRasterBandH *huge_msi_img, // satellite msi image
        int pd, int pdm                // dim1 = enz et dim2 = msi 
        );

double norm_with_t_threshold(
        double *v_sun,       // theoric sun information (seen ? utm : NAN)
        double *out_pan,     // vertices intensity
        int nv,              // number of vertices
        double t             // threshold
        );

double get_shadow_threshold(
        double *v_sun,               // theoric sun information (seen ? utm : NAN)
        double *out_pan,             // vertices intensity
        int nv,                      // number of vertices
        double t_min, double t_max   // threshold bound
        );

void create_shadow_mask(
        char *filename_sun,          // file with theoric sun information
        double *out_pan,             // vertices intensity
        double *out_sun,             // vertices lighting (shadow ? 0 : 1)
        double t_min, double t_max   // threshold bound
        );

#endif //COLORIZE_VERTICES_FROM_ONE_IMAGE_H
