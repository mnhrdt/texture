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
#include "pickopt.h"


int main(int c, char *v[])
{
    double ox = atof(pick_option(&c, &v, "ox", "000"));
    double oy = atof(pick_option(&c, &v, "oy", "000"));
    double oz = atof(pick_option(&c, &v, "oz", "000"));
    int line = atof(pick_option(&c, &v, "l", "0"));

    if (c < 4)
        return fprintf(stderr, "usage:\n\t"
                "%s mesh.off vc.png out.ply\n",*v);
                //0 1        2      3           
                
    char *filename_mesh = v[1];
    char *filename_colors = v[2];
    char *filename_ply = v[3];

    struct trimesh m;
    trimesh_read_from_off(&m, filename_mesh);

    int w, h, pd;
    double *vc = iio_read_image_double_vec(filename_colors, &w, &h, &pd);

    double origin[3] = {ox, oy, oz};
    trimesh_write_to_coloured_ply(filename_ply, &m, vc+pd*line*w, origin); 

    return 0;
}

