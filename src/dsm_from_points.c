#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "trimesh.c"
#include "iio.h"

// TODO : problème si resolution négative

void get_dsm_size(struct trimesh *pc, double res[2], double wm[2], double hm[2])
{
    double whm[4] = {INFINITY, -INFINITY, INFINITY, -INFINITY};

    for (int i = 0; i < pc->nv; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            whm[2*j] = (whm[2*j] > pc->v[3*i+j]/res[j]) ? pc->v[3*i+j]/res[j] : whm[2*j];
            whm[2*j+1] = (whm[2*j+1] < pc->v[3*i+j]/res[j]) ? pc->v[3*i+j]/res[j] : whm[2*j+1];
        }
    }
    printf("x %f %f y %f %f\n", whm[0],whm[1],whm[2],whm[3]);
    for (int j = 0; j < 2; j++)
    {
        wm[j] = floor(whm[j]);
        hm[j] = floor(whm[j+2]);
    }
}

void fill_dsm_histo(struct trimesh *pc, double res[2], double wm[2], double hm[2], double *dsm_histo)
{
    int w = (int) (wm[1] - wm[0]) +1;
    int h = (int) (hm[1] - hm[0]) +1;
    for (int nv = 0; nv < pc->nv; nv++)
    {
        int i = (int) (floor(pc->v[3*nv]/res[0]) - wm[0]);
        int j = (int) (floor(pc->v[3*nv+1]/res[1]) - hm[0]);
        dsm_histo[2*(i+j*w)] += pc->v[3*nv+2];
        dsm_histo[2*(i+j*w)+1] += 1; 
    }
    for (int i = 0; i < w*h; i++)
        dsm_histo[2*i] /= fmax(dsm_histo[2*i+1], 1.0);
}

int main(int c, char *v[])
{
    if (c < 4)
        return fprintf(stderr, "usage: \n\t"
                "%s cloud.off res_x res_y out.tif\n", *v);
    char *filename_cloud = v[1];
    double res[2] = {atof(v[2]), atof(v[3])};
    char *filename_dsm = v[4];
    
    struct trimesh pc;
    trimesh_read_from_off(&pc, filename_cloud);
    
    double wm[2] = {INFINITY, -INFINITY};
    double hm[2] = {INFINITY, -INFINITY};
    get_dsm_size(&pc, res, wm, hm);
    int w = (int) (wm[1] - wm[0]) +1;
    int h = (int) (hm[1] - hm[0]) +1;
    printf("w %d h %d\n wm %f %f hm %f %f\n", w, h, wm[0], wm[1], hm[0], hm[1]);
    if (w < 1 || h < 1)
        return fprintf(stderr, "ERROR: dsm size nonpositive, w = %d, h = %d\n", w, h);

    double *dsm_histo = malloc(2*w*h*sizeof(double));

    fill_dsm_histo(&pc, res, wm, hm, dsm_histo);

    double *dsm = malloc(w*h*sizeof(double));
    for (int i = 0; i < w*h; i++)
        dsm[i] = dsm_histo[2*i];

    iio_save_image_double(filename_dsm, dsm, w, h);
}
