#include <assert.h>
#include <stdio.h>
#include <gdal.h>
#include <cpl_conv.h>
#include <libgen.h>
#include <math.h>
#include <string.h>

#include "rpc.c"
#include "fail.c"
#include "iio.h"
#include "trimesh.c"

int utm_from_lonlat(double out_eastnorth[2], double lon, double lat);
void lonlat_from_eastnorthzone(double out_lonlat[2], double e, double n, int z);

// extrapolate by nearest value
static float getpixel_f(float *I, int w, int h, int i, int j)
{
    if (i < 0) i = 0;
    if (j < 0) j = 0;
    if (i >= w) i = w-1;
    if (j >= h) j = h-1;
    return I[i+j*w];
}

// some useful functions
double scalar_product(double *u, double *v, int dim)
{
    double sp = 0;
    for (int i = 0; i < dim; i++)
        sp += u[i]*v[i];
    return sp;
}

void cross_product(double axb[3], double a[3], double b[3])
{
    // a0 a1 a2
    // b0 b1 b2
    axb[0] = a[1]*b[2] - a[2]*b[1];
    axb[1] = a[2]*b[0] - a[0]*b[2];
    axb[2] = a[0]*b[1] - a[1]*b[0];
}

double euclidean_norm(double *x, int n)
{
    return n > 0 ? hypot(*x, euclidean_norm(x+1, n-1)) : 0;
}

static bool is_in_crop_int(int *ij, int *xywh)
{
    return ((ij[0]>= xywh[0]) && (ij[0] < xywh[0] + xywh[2]) 
            && (ij[1]>= xywh[1]) && (ij[1] < xywh[1] + xywh[3])); 
}
// conversion from hsv to rgb (h in [0, 360], s v r g b in [0, 1])
void hsv2rgb(double hsv[3], double rgb[3])
{
    int hi = (int) floor(hsv[0]/60) % 6;
    double f = hsv[0]/60 - hi;
    double l = hsv[2] * (1 - hsv[1]);
    double m = hsv[2] * (1 - f * hsv[1]);
    double n = hsv[2] * (1 - (1 - f) * hsv[1]);
    if (hi == 0)
    {
        rgb[0] = hsv[2];
        rgb[1] = n;
        rgb[2] = l;
    }
    if (hi == 1)
    {
        rgb[0] = m;
        rgb[1] = hsv[2];
        rgb[2] = l;
    }
    if (hi == 2)
    {
        rgb[0] = l;
        rgb[1] = hsv[2];
        rgb[2] = n;
    }
    if (hi == 3)
    {
        rgb[0] = l;
        rgb[1] = m;
        rgb[2] = hsv[2];
    }
    if (hi == 4)
    {
        rgb[0] = n;
        rgb[1] = l;
        rgb[2] = hsv[2];
    }
    if (hi == 5)
    {
        rgb[0] = hsv[2];
        rgb[1] = l;
        rgb[2] = m;
    }
}

void create_colormap(char *filename_map, int nimages, int w, int h) // nimages max = 48
{
    float *map = malloc(6 * w * h * nimages * sizeof(float));
    double Hsv[8] = {0, 45, 60, 120, 180, 225, 270, 315};
    double hSv[3] = {1, 0.33, 0.67};
    double hsV[2] = {1, 0.7};
    for (int ni = 0; ni < nimages; ni++)
    {
        double hsv[3] = {Hsv[ni%8], hSv[(ni/16)%3], hsV[(ni/8)%2]};
        double rgb[3];
        hsv2rgb(hsv, rgb);
        for (int i = ni*w; i < (ni+1)*w; i++)
            for (int j = 0; j < h; j++)
                for (int l = 0; l < 3; l++)
                    map[3*(i+j*nimages*w)+l] = 255*rgb[l];
    }
    for (int i = 3*w*h*nimages; i < 6*w*h*nimages; i++)
        map[i] = 0;
    iio_save_image_float_vec(filename_map, map, nimages*w, 2*h, 3);
    free(map);

}

// coordinates of the normal to a triangle in a 3D space
void triangle_normal(double n[3], double a[3], double b[3], double c[3]) // les sommets sont donnés dans le sens direct
{
    double u[3];
    double v[3];
    for (int i = 0; i < 3; i++) u[i] = b[i] - a[i];
    for (int i = 0; i < 3; i++) v[i] = c[i] - a[i];
    cross_product(n, u, v);

    double norm = euclidean_norm(n, 3);
    for (int i = 0; i < 3; i++) n[i] /= norm;

    norm = euclidean_norm(n, 3);
    if (norm < 0.9999 || norm > 1.0000001)
        printf("WARNING: normalisation error in triangle_normal, norme = %.16lf\n", norm);
}

// check if (xmin < x < xmax) && (ymin < y < ymax)
bool xy_are_in_bounds(double *xy, // tableau de coords (x1, x2,..., y1,y2,..)
        int l, // longueur des x/y
        double xmin, double xmax, double ymin, double ymax)
{       
    bool a = true;
    for (int i = 0; i < l; i++) 
        a &= (xy[i] >= xmin && xy[i] <= xmax &&
                xy[i+l] >= ymin && xy[i+l] <=ymax);
    return a;
}



// get the satellite direction using only the rpc data
void camera_direction(double n[3], struct rpc *r)
{
    // initialise height and fill 3rd vector coordinate
    double z = 0;
    n[2] = -1;

    // get first 3D point using localisation 
    double ijh[3] = {500, 500, z};
    double lonlat[2] = {0, 0};
    rpc_localization(lonlat, r, ijh);
    double en[2];
    utm_from_lonlat(en, lonlat[0], lonlat[1]);
    for (int i = 0; i < 2; i++)
        n[i] = en[i];

    // get second 3D point combining localisation and projection
    double ij[2] = {0, 0};
    double lonlatheight[3] = {lonlat[0], lonlat[1], z};
    rpc_projection(ij, r, lonlatheight);
    for (int i = 0; i < 2; i++)
        ijh[i] = ij[i];
    ijh[2] = z + 1;
    rpc_localization(lonlat, r, ijh);
    utm_from_lonlat(en, lonlat[0], lonlat[1]);

    // fill in the first two vector coordinates
    for (int i = 0; i < 2; i++)
        n[i] -= en[i];

    // normalise direction vector
    double norm = euclidean_norm(n, 3);
    for (int i = 0; i < 3; i++) n[i] /= norm;

    norm = euclidean_norm(n, 3);
    if (norm < 0.9999 || norm > 1.0000001)
        printf("WARNING: normalisation error in camera_direction, norme = %.16lf\n", norm);
}


// elevate a georeferenced DSM to an ascii point cloud

// elevate a georeferenced DSM to an ascii point cloud
// optionally, colorize the point cloud from a given reference image
// optionally, create a ply file with the desired connectivity
#include "pickopt.c"

// main_colorize: like elevate, but produce a color image
//
// input_1 : a georeferenced DSM of size WxH
// input_2 : a georeferenced d-band image (of any size, typically huge)
// input_3 : the RPC of the d-band image
// output  : a d-band image of size WxH, with colors from the correct place
//
// note: formally, for the input we have to give the signed zone
//
int main_colormultiple(int c, char *v[])
{
    double ox = atof(pick_option(&c, &v, "ox", "0"));
    double oy = atof(pick_option(&c, &v, "oy", "0"));
    float offset[2] = {ox, oy};
    if (c < 5)
        return fprintf(stderr, "usage:\n\t"
                "%s mesh.ply img_i.tif match_i.tif out.ply\n",*v);
    //0 1        2         3           4      
    char *filename_mesh = v[1];
    char *filename_im = v[2];
    char *filename_m = v[3];
    char *filename_ply = v[4];

    printf("début\n");
    int wi, hi;
    float *im = iio_read_image_float(filename_im, &wi, &hi);
    if (!im)
        return fprintf(stderr, "iio_read(%s) failed\n", filename_im);
    printf("image chargée\n");

    int wh, un, pd;
    double *match = iio_read_image_double_vec(filename_m, &wh, &un, &pd);
    if (!match)
        return fprintf(stderr, "iio_read(%s) failed\n", filename_m);
    printf("matchs chargés\n");

    struct trimesh m;
    trimesh_read_from_off(&m, filename_mesh);
    printf("trimesh chargé\n");


    if (m.nv != wh)
        return fprintf(stderr, "mesh and matches dimensions mismatch : m.nv %d wh %d \n", m.nv, wh);

    double *vc = malloc(3 * m.nv * sizeof(double));
    printf("création vecteur couleur\n");

    for (int nv = 0; nv < m.nv; nv++)
    {
        int ij[2];
        int xywh[4] = {0, 0, wi, hi};
        for (int i = 0; i < 2; i++)
            ij[i] = (int) round(match[2*nv+i] + offset[i]);
        if (is_in_crop_int(ij, xywh))
            for (int k = 0; k < 3; k++)
                vc[3*nv+k] = im[ij[0]+(ij[1])*wi];
        else 
        {
            vc[3*nv+0] = NAN;
            vc[3*nv+1] = NAN;
            vc[3*nv+2] = NAN;
        }
    }

    printf("avant écriture\n");
    trimesh_write_to_coloured_ply(filename_ply, &m, vc); 

    printf("fini !\n");
    return 0;
}

int main(int c, char *v[]) { return main_colormultiple(c,v); }

