#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "iio.h"

// extrapolate by nearest neighbour
static void keep_in_bounds(int w, int h, int *ij)
{
                if (ij[0] < 0) ij[0] = 0;
                if (ij[1] < 0) ij[1] = 0;
                if (ij[0] >= w) ij[0] = w-1;
                if (ij[1] >= h) ij[1] = h-1;
}
    
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


int main_colorize_with_shadows(int c, char *v[])
{
        char *filename_dsm = v[1];
        char *filename_img = v[2];
        char *filename_P = v[3];
        char *filename_out = v[4];

        // read the whole input DSM (typically, rather small)
        int w, h;
        float *x = iio_read_image_float(filename_dsm, &w, &h);
        if (!x)
                return fprintf(stderr, "iio_read(%s) failed\n", filename_dsm);

        // read the input image (small jpg, png or tif crop corresponding to the DSM)
        int wi, hi;
        float *img = iio_read_image_float(filename_img, &wi, &hi);
        if (!img)
                return fprintf(stderr, "iio_read(%s) failed\n", filename_img);

        // allocate space of same size and initialize to -100
        float *img_copy = malloc(wi * hi * sizeof(float));
        for (int i=0; i < wi*hi; i++)
                img_copy[i] = -100;

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

        // loop over lidar points and fill img_copy with height of lidar point
        // keep only the highest heigtht for each point to deal with occlusions
        for (int j = 0; j < h; j++)
                for (int i = 0; i < w; i++)
                {
                        double z = x[j*w+i];
                        if (z<0)
                                z = 20.0;
                        double xyz1[4] = {i, j, z, 1};
                        double ij_approx[3];
                        
                        matrix_product_mxn_nxp(ij_approx, P, xyz1, 3, 4, 1); 

                        int ij_int[2];
                        for (int k = 0; k < 2; k++)
                                ij_int[k] = (int) round(ij_approx[k]);

                        keep_in_bounds(wi, hi, ij_int);
                        if (img_copy[ij_int[1]*wi+ij_int[0]] < z)
                                img_copy[ij_int[1]*wi+ij_int[0]] = z;
                }
        

        // allocate space for output
        float *out_x = malloc(2 * w * h * sizeof(float));

        // loop over lidar point and give them coordinate of projection on image
        // only if there are no occlusions. Otherwise nan.
        for (int i = 0; i < w; i++)
                for (int j = 0; j < h; j++)
                {
                        double z = x[j * w+i];
                        if (z<0)
                                z = 20.0;
                        double xyz1[4] = {i, j, z, 1};
                        double ij_approx[3];
                        
                        matrix_product_mxn_nxp(ij_approx, P, xyz1, 3, 4, 1);

                        int ij_int[2];
                        for (int k = 0; k < 2; k++)
                                ij_int[k] = (int) round(ij_approx[k]);

                        keep_in_bounds(wi, hi, ij_int);
                        if (img_copy[ij_int[1]*wi + ij_int[0]] == z)
                                for (int k = 0; k < 2; k++)
                                        out_x[2*(j*w+i)+k] = ij_approx[k];
                        else
                                for (int k = 0; k < 2; k++)
                                        out_x[2*(j*w+i)+k] = NAN;
                }

        iio_save_image_float_vec(filename_out, out_x, w, h, 2);

        return 0;
}

int main(int c, char *v[]) { return main_colorize_with_shadows(c,v); }
