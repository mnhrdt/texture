#include <assert.h>
#include <stdio.h>
#include <gdal.h>
#include <cpl_conv.h>
#include <libgen.h>
#include <math.h>
#include <string.h>

#include "rpc.h"
#include "fail.c"
#include "iio.h"
#include "trimesh.h"
#include "normals.h"
#include "colorize_vertices_from_one_image.h"
#include "pickopt.h"

void lonlat_from_eastnorthzone(
        double out_lonlat[2],            // longitude and latitude
        double e,                        // easting 
        double n,                        // northing
        int z)                           // signed zone
    ;

// cubic interpolation in dimension 1
static float cubic_interpolation(float v[4], float x)
{
        return v[1] + 0.5 * x*(v[2] - v[0]
                        + x*(2.0*v[0] - 5.0*v[1] + 4.0*v[2] - v[3]
                        + x*(3.0*(v[1] - v[2]) + v[3] - v[0])));
}

// bicubic interpolation in dimension 2 (needs a cell of size 4x4)
static float bicubic_interpolation_cell(float *p, float x, float y)
{
        float v[4];
        v[0] = cubic_interpolation(p + 4*0, y);
        v[1] = cubic_interpolation(p + 4*1, y);
        v[2] = cubic_interpolation(p + 4*2, y);
        v[3] = cubic_interpolation(p + 4*3, y);
        return cubic_interpolation(v, x);
}

// eval a gdal raster using bicubic interpolation
static float gdal_getpixel_bicubic(GDALRasterBandH img, double x, double y)
{
        // caching of a silly malloc
        static float *roi = NULL;
        if (!roi) roi = CPLMalloc(4*4*sizeof*roi);

        // TODO: some sort of bilinear or bicubic interpolation ?
        x -= 1;
        y -= 1;
        int i = floor(x);
        int j = floor(y);
        int r = GDALRasterIO(img, GF_Read, i,j,4, 4, roi,4,4, GDT_Float32, 0,0);
        (void) r;
        return bicubic_interpolation_cell(roi, y - j, x - i);
//        return roi[5];
}

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
        int pd, int pdm)               // dim1 = enz et dim2 = msi 
{ 
    printf("pdm = %d\n", pdm);
    // initialize outputs to values if not seen by camera
    for (int i = 0; i < 3 * m->nv; i++){
        out_rgb[i] = NAN;
        out_pan[i] = NAN;}
    for (int i = 0; i < pdm * m->nv; i++)
        out_msi[i] = NAN;

    // for each vertex, get color informations if pixel seen by satellite
    for (int v = 0; v < m->nv; v++)
        if (!isnan(coord[pd*v]) && coord[pd*v+2] > -500)
        {
            // get lonlat coordinates from utm coordinates
            double lonlat[2] = {0, 0};
            lonlat_from_eastnorthzone(lonlat, coord[pd*v], coord[pd*v+1], signed_zone);
            double lonlatheight[3] = {lonlat[0], lonlat[1], coord[pd*v+2]};

            // get projection coordinates in pan and msi satellite images
            double ij_pan[2];
            double ij_msi[2];
            rpc_projection(ij_pan, pan_rpc, lonlatheight);
            rpc_projection(ij_msi, msi_rpc, lonlatheight);

            // fill panchromatic output
            double intensity = gdal_getpixel_bicubic(huge_pan_img[0], ij_pan[0], ij_pan[1]);
            for (int i = 0; i < 3; i++)
                out_pan[3*v+i] = intensity;

            // fill multispectral output
            double msi[pdm];
            for (int l = 0; l < pdm; l++){
                msi[l] = gdal_getpixel_bicubic(huge_msi_img[l], ij_msi[0], ij_msi[1]);
                if (isnan(msi[l]))
                    printf("vertex with nan %d\n", v);
                out_msi[pdm * v + l] = msi[l];}
            if (v == 955649)
                printf("msi %lf %lf %lf %lf %lf %lf %lf %lf \n", out_msi[pdm * v + 0], out_msi[pdm * v + 1], out_msi[pdm * v + 2], out_msi[pdm * v + 3], out_msi[pdm * v + 4], out_msi[pdm * v + 5], out_msi[pdm * v + 6], out_msi[pdm * v + 7]);

            // fill rgb output obtained by pansharpening
            double rgb[3] = {msi[4], 0.8 * msi[2] + 0.1 * msi[5], 1.2 * msi[1]};
            double a = intensity/(rgb[0] + rgb[1] + rgb[2]);
            for (int i = 0; i < 3; i++){
                out_rgb[3 * v + i] = a * rgb[i];
                if (out_rgb[3 * v + i] > 10000)
                    printf("vertex with rgb problem %d\n", v);}
            if (v == 955649)
                printf("rgb %lf %lf %lf pan %lf \n", out_rgb[3 * v + 0], out_rgb[3 * v + 1], out_rgb[3 * v + 2], a);
        }
}

double norm_with_t_threshold(
        double *v_sun,       // theoric sun information (seen ? utm : NAN)
        double *out_pan,     // vertices intensity
        int nv,              // number of vertices
        double t)            // threshold
{
    // initialize number of differences between theory and reality
    double norm = 0;

    // theoric information coherent with thresholding, norm +=1
    for (int v = 0; v < nv; v++)
    {
        if (isnan(out_pan[3*v]))      // check if vertex seen by camera
                continue;
        if (((out_pan[3*v] > t) && isnan(v_sun[3 * v])) || ((out_pan[3*v] <= t) && !isnan(v_sun[3*v])))
            norm += 1;
    }
    return norm;
}

double get_shadow_threshold(
        double *v_sun,               // theoric sun information (seen ? utm : NAN)
        double *out_pan,             // vertices intensity
        int nv,                      // number of vertices
        double t_min, double t_max)  // threshold bound
{
    // initialize threshold and l1 norms with this threshold
    double t = t_min;
    double l1norm = norm_with_t_threshold(v_sun, out_pan, nv, t);
    double l1norm_old = l1norm;

    // keep incrementing while the norm decreases
    while (l1norm_old >= l1norm && t < t_max)
    {
        l1norm_old = l1norm;
        t += 1;
        l1norm = norm_with_t_threshold(v_sun, out_pan, nv, t);
    }
    return t;
}

void create_shadow_mask(
        char *filename_sun,          // file with theoric sun information
        double *out_pan,             // vertices intensity
        double *out_sun,             // vertices lighting (shadow ? 0 : 1)
        double t_min, double t_max)  // threshold bound
{
    // read theoric sun information (seen by sun ? utm : NAN)
    int nv, un, pd;
    double *v_sun = iio_read_image_double_vec(filename_sun, &nv, &un, &pd);

    // get shadow threshold and fill practical shadow location
    double t = get_shadow_threshold(v_sun, out_pan, nv, t_min, t_max);
    printf("threshold t = %lf", t);
    for (int v = 0; v < 3 * nv; v++)
        out_sun[v] = (out_pan[v] > t) ? 1 : 0;
}

int main(int c, char *v[])
{
    double t_min = atof(pick_option(&c, &v, "tmin", "0"));
    double t_max = atof(pick_option(&c, &v, "tmax", "1000"));
    if (c < 5)
        return fprintf(stderr, "usage:\n\t"
                "%s mesh.off pan.ntf pan.xml msi.ntf msi.xml"
                "utm_coord.tif vc.tif\n", *v);

    char *filename_mesh = v[1];       // scaled mesh
    char *filename_pan = v[2];        // panchromatic satellite image .ntf
    char *filename_pan_rpc = v[3];    // pan sat image rpc
    char *filename_msi = v[4];        // multispectral satellite image .ntf
    char *filename_msi_rpc = v[5];    // msi sat image rpc
    int signed_zone = atoi(v[6]);     // utm signed zone
    char *filename_coord = v[7];      // vertices utm coords if seen by sat
    char *filename_sun = v[8];        // vertices utm coords if seen by sun
    char *filename_out_pan = v[9];    // vertices intensity
    char *filename_out_msi = v[10];    // vertices msi values
    char *filename_out_rgb = v[11];   // vertices rgb from pansharpening
    char *filename_out_sun = v[12];   // vertices shadow location

    GDALAllRegister();

    // open the reference image and obtain its pixel dimension 
    GDALDatasetH huge_dataset = GDALOpen(filename_pan, GA_ReadOnly);
    int pdm = GDALGetRasterCount(huge_dataset);
    GDALRasterBandH huge_pan_img[pdm];
    for (int i = 0; i < pdm; i++)
        huge_pan_img[i] = GDALGetRasterBand(huge_dataset, i+1);

    // open the reference rpc
    struct rpc pan_rpc[1];
    read_rpc_file_xml(pan_rpc, filename_pan_rpc);

    // open the reference image and obtain its pixel dimension 
    huge_dataset = GDALOpen(filename_msi, GA_ReadOnly);
    pdm = GDALGetRasterCount(huge_dataset);
    GDALRasterBandH huge_msi_img[pdm];
    for (int i = 0; i < pdm; i++)
        huge_msi_img[i] = GDALGetRasterBand(huge_dataset, i+1);

    // open the reference rpc
    struct rpc msi_rpc[1];
    read_rpc_file_xml(msi_rpc, filename_msi_rpc);

    // read vertices utm coordinates
    int nv, un, pd;
    double *utm_coord = iio_read_image_double_vec(filename_coord, &nv, &un, &pd);
    if (!utm_coord)
        return fprintf(stderr, "iio_read(%s) failed\n", filename_coord);

    // read scaled mesh
    struct trimesh m[1];
    trimesh_read_from_off(m, filename_mesh);

    // create vertices output informations
    double *out_pan = malloc(3 * m->nv * sizeof(double));
    double *out_sun = malloc(3 * m->nv * sizeof(double));
    double *out_rgb = malloc(3 * m->nv * sizeof(double));
    double *out_msi = malloc(pdm * m->nv * sizeof(double));

    // get intensity, msi and pansharpened rgb values
    colorize(out_pan, out_msi, out_rgb, m, utm_coord, pan_rpc, msi_rpc, signed_zone, huge_pan_img, huge_msi_img, pd, pdm);

    // get shadow location
    create_shadow_mask(filename_sun, out_pan, out_sun, t_min, t_max);

    // save outputs
    iio_save_image_double_vec(filename_out_pan, out_pan, nv, 1, 3);
    iio_save_image_double_vec(filename_out_sun, out_sun, nv, 1, 3);
    iio_save_image_double_vec(filename_out_rgb, out_rgb, nv, 1, 3);
    iio_save_image_double_vec(filename_out_msi, out_msi, nv, 1, pdm);

    // free allocated memory
    free(out_pan); free(out_msi); free(out_rgb); free(out_sun);

    return 0;
}


