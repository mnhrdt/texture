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
    double ox = atof(pick_options(&c, &v, "ox", "0"));
    double oy = atof(pick_options(&c, &v, "oy", "0"));
    double oz = atof(pick_options(&c, &v, "oz", "0"));
    if (c < 3)
        return fprintf(stderr, "usage:\n\t"
                "%s dsm.tif zone mesh.off\n", *v);
                //0 1       2    3    
    char *filename_dsm = v[1];
    int signed_zone    = atoi(v[2]);
    char *filename_mesh = v[3];

    // read the whole input DSM (typically, rather small)
    int w, h;
    float *x = iio_read_image_float(filename_dsm, &w, &h);
    if (!x)
        return fprintf(stderr, "iio_read(%s) failed\n", filename_dsm);

    // read georeferencing transform using GDAL
    double origin[2] = {0, 0};
    double scale[2] = {1, 1};

    GDALAllRegister();
    GDALDatasetH gdal_dataset = GDALOpen(filename_dsm, GA_ReadOnly);
    if (gdal_dataset == NULL)
        fprintf(stderr, "GDALOpen(%s) failed\n", filename_dsm);
    double tmp[6];
    if (GDALGetGeoTransform(gdal_dataset, tmp) == CE_None) {
        origin[0] = tmp[0], origin[1] = tmp[3];
        scale[0] = tmp[1], scale[1] = tmp[5];
    } else {
        fprintf(stderr, "WARNING: not found origin and scale info\n");
    }

    struct trimesh m[1];
    trimesh_create_from_dem_with_offset(m, x, w, h, ox/0.3, oy/0.3, oz);

    trimesh_write_to_off(filename_mesh, m);

    return 0;
}

