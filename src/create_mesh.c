#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gdal.h>

#include "iio.h"
#include "trimesh.h"
#include "rpc.h"

int main(int c, char *v[])
{
    if (c < 4)
        return fprintf(stderr, "usage:\n\t"
                "%s dsm.tif mesh.off mesh.ply\n", *v);
                //0 1       2    3    
    char *filename_dsm = v[1];
    char *filename_off = v[2];
    char *filename_ply = v[3];

    // read the whole input DSM (typically, rather small)
    int w, h;
    float *x = iio_read_image_float(filename_dsm, &w, &h);
    if (!x)
        return fprintf(stderr, "iio_read(%s) failed\n", filename_dsm);

     // read georeferencing transform using GDAL
    double scale[3] = {1, 1, 1};
    GDALAllRegister();
    GDALDatasetH gdal_dataset = GDALOpen(filename_dsm, GA_ReadOnly);
    if (gdal_dataset == NULL)
        fprintf(stderr, "GDALOpen(%s) failed\n", filename_dsm);
    double tmp[6];
    if (GDALGetGeoTransform(gdal_dataset, tmp) == CE_None) {
        scale[0] = tmp[1], scale[1] = tmp[5];
    } else {
        fprintf(stderr, "WARNING: not found origin and scale info\n");
    }

    struct trimesh m[1];
    trimesh_create_from_dem_with_scale(m, x, w, h, scale);

    trimesh_write_to_ply(filename_ply, m);
    trimesh_write_to_off(filename_off, m);

    return 0;
}

