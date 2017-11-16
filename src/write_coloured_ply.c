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
/*void mean(double k)
{
    k  = 0;
}
void vc_mean2(double *vcf, double *vcs, int nv, int nb_vc)
{
    int count;
    double vc[3*nb_vc]; 
    for (int j = 0; j < nv; j++)
    {
        count = 0;
        for (int i = 0; i < 3*nb_vc; i++)
            vc[i] = NAN;
        for (int i = 0; i < nb_vc; i++)
            if (!isnan(vcs[3*(i*nv + j)]))
            {
                for (int k = 0; k < 3; k++)
                    vc[3*count + k] = vcs[3*(i*nv+j)+k];
                count += 1;
            }
        for (int k = 0; k < 3; k++)
            vcf[3*j +k] = sum[k]/count;
    }
}
*/
void vc_mean(double *vcf, double *vcs, int nv, int nb_vc)
{
    double sum[3];
    int count;
    for (int j = 0; j < nv; j++)
    {
        count = 0;
        for (int k = 0; k < 3; k++)
            sum[k] = 0;
        for (int i = 0; i < nb_vc; i++)
            if (!isnan(vcs[3*(i*nv + j)]))
            {
                count +=1;
                for (int k = 0; k < 3; k++)
                    sum[k] += vcs[3*(i*nv + j) + k];
            }
        for (int k = 0; k < 3; k++)
            vcf[3*j +k] = sum[k]/count;
    }
}


int main(int c, char *v[])
{
    double t = atof(pick_option(&c, &v, "t", "255"));

    if (c < 4)
        return fprintf(stderr, "usage:\n\t"
                "%s mesh.off vc_i.tif out.ply\n",*v);
                //0 1        2         3           
    printf("c %d\n", c);
    int nb_vc = c-3;
    char *filename_mesh = v[1];
    char *filename_ply = v[2 + nb_vc];

    struct trimesh m;
    trimesh_read_from_off(&m, filename_mesh);
    printf("trimesh chargé\n");

    double *vcs;
    vcs = malloc(3 * m.nv * nb_vc * sizeof(double));
    for (int i = 0; i < nb_vc; i++)
    {
        int nv, one, pd;
        double *vc = iio_read_image_double_vec(v[2+i], &nv, &one, &pd);
        if (!vc || one != 1 || pd != 3)
            return fprintf(stderr, "iio_read(%s) failed\n", v[2+i]);
        if (m.nv != nv)
            return fprintf(stderr, "mesh and matches dimensions mismatch :" 
                    "m.nv %d w %d \n", m.nv, nv);
        for (int j = 0; j < 3*nv; j++)
            vcs[j+i*3*nv] = vc[j]; 
    }

    double *vcf;
    vcf = malloc(3 * m.nv * sizeof(double));

    if (nb_vc == 1)
        for (int i = 0; i < 3 * m.nv; i++)
            vcf[i] = vcs[i];
    else
        vc_mean(vcf, vcs, m.nv, nb_vc);

    printf("avant écriture\n");
    double origin[3] = {0, 0, 0};
    trimesh_write_to_coloured_ply2(filename_ply, &m, vcf, origin); 

    free(vcf); free(vcs);
    printf("fini !\n");
    return 0;
}


