#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gdal.h>

#include "iio.h"
#include "drawtriangle.c"

static void matrix_product_mxn_nxp(double *ab, double *a, double *b, int m, 
                int n, int p)
{
        for (int i = 0; i < m; i++)
                for (int j = 0; j < p; j++)
                {
                        ab[i*p+j] = 0;
                        for (int k = 0; k <n; k++)
                                ab[i*p+j] += a[i*n+k] * b[k*p+j];
                }
}


int main_vector_proj(int c, char *v[])
{
        char *filename_P = v[1];
        char *filename_vec = v[2];

        // read P projection matrix of lidar on cropped image
        FILE *proj;
        proj = fopen(filename_P,"r");
        if (!proj)
                return fprintf(stderr, "fopen(%s) failed\n", filename_P);

        double P[12];

        for (int i = 0; i < 12; i++)
        {
                int r = fscanf(proj, "%lf", &P[i]);
		if (r!=1)
			return fprintf(stderr, "could not read element %d of P", i);
        }
        fclose(proj);

        FILE *bias;
        bias = fopen(filename_vec,"r");
        if (!bias)
                return fprintf(stderr, "fopen(%s) failed\n", filename_vec);
        double xysz[4];
        for (int i = 0; i < 4; i++)
        {
                int r = fscanf(bias, "%lf", &xysz[i]);
		if (r!=1)
			return fprintf(stderr, "could not read element %d of %s", i, filename_vec);
        }
        
        double xyz1[4] = {xysz[0], xysz[1], xysz[3], 0};

        double ij1[3] = {0};
        
        matrix_product_mxn_nxp(ij1, P, xyz1, 3, 4, 1);

        printf("%lf %lf\n", ij1[0], ij1[1]);


        return 0;
}

int main(int c, char *v[]) { return main_vector_proj(c,v); }
