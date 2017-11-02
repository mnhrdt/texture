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
#include "pickopt.c"

int main(int c, char *v[])
{
    double t = atof(pick_option(&c, &v, "t", "800"));

    if (c < 4)
        return fprintf(stderr, "usage:\n\t"
                "%s mesh.off vc.tif out.ply\n",*v);
                //0 1        2      3            
    char *filename_mesh = v[1];
    char *filename_vc = v[2];
    char *filename_ply = v[3];

    int nv, one, pd;
    double *vc = iio_read_image_double_vec(filename_vc, &nv, &one, &pd);
    if (!vc || one != 1 || pd != 3)
        return fprintf(stderr, "iio_read(%s) failed\n", filename_vc);
    printf("vecteur de couleurs chargé\n");

    struct trimesh m;
    trimesh_read_from_off(&m, filename_mesh);
    printf("trimesh chargé\n");
    if (m.nv != nv)
        return fprintf(stderr, "mesh and matches dimensions mismatch :" 
                "m.nv %d w %d \n", m.nv, nv);

    printf("avant écriture\n");
    trimesh_write_to_coloured_ply(filename_ply, &m, vc, t); 

    printf("fini !\n");
    return 0;
}


