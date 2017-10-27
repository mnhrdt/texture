#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gdal.h>

#include "iio.h"
#include "trimesh.c"
#include "rpc.c"
#include "pickopt.c"

int main(int c, char *v[])
{
    int ox = atof(pick_option(&c, &v, "ox", "0"));
    int oy = atof(pick_option(&c, &v, "oy", "0"));
    double oz = atof(pick_option(&c, &v, "oz", "0"));
    if (c < 4)
        return fprintf(stderr, "usage:\n\t"
                "%s dsm.tif zone mesh.off mesh.ply\n", *v);
                //0 1       2    3    
    char *filename_dsm = v[1];
    char *filename_off = v[2];
    char *filename_ply = v[3];

    // read the whole input DSM (typically, rather small)
    int w, h;
    float *x = iio_read_image_float(filename_dsm, &w, &h);
    if (!x)
        return fprintf(stderr, "iio_read(%s) failed\n", filename_dsm);

    struct trimesh m[1];
    trimesh_create_from_dem_with_offset(m, x, w, h, ox, oy, oz);

    trimesh_write_to_ply(filename_ply, m);
    trimesh_write_to_off(filename_off, m);

    return 0;
}

