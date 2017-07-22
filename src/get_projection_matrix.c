#include <assert.h>
#include <stdio.h>
#include <gdal.h>
#include <cpl_conv.h>

#include "rpc.c"
#include "fail.c"
#include "iio.h"

int utm_from_lonlat(double out_eastnorth[2], double lon, double lat);
void lonlat_from_eastnorthzone(double out_lonlat[2], double e, double n, int z);

static float gdal_getpixel(GDALRasterBandH img, double pi, double pj)
{
	// caching of a silly malloc
	static float *roi = NULL;
	if (!roi) roi = CPLMalloc(1*1*sizeof*roi);

	// TODO: some sort of bilinear or bicubic interpolation ?
	int i = round(pi);
	int j = round(pj);
	int r = GDALRasterIO(img, GF_Read, i,j,1,1, roi,1,1, GDT_Float32, 0,0);
	return roi[0*0+0];
}


void eval_proj(double *ij, double e, double n, double z,  
                struct rpc *huge_rpc, int signed_zone)
{
        // compute lonlat from eastnorth
        double lonlat[3] = {0, 0, z};
        lonlat_from_eastnorthzone(lonlat, e, n, signed_zone);

        // compute coordiantes in huge image
        rpc_projection(ij, huge_rpc, lonlat);
}

static void matrix_product_mxn_nxp(double *ab, double *a, double *b, 
                int m, int n, int p)
{
        // a, b and ab filled with juxtaposed lines, a : m lines and n columns
        for (int i = 0; i < m; i++)
                for (int j = 0; j < p; j++)
                {
                        ab[i*p+j] = 0;
                        for (int k = 0; k < n; k++)
                                ab[i*p+j] += a[i*n+k] * b[k*p+j];
                }
}
// elevate a georeferenced DSM to an ascii point cloud
//
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
int main_get_projection_matrix(int c, char *v[])
{
	double offset_x = atof(pick_option(&c, &v, "-offset_x", "0"));
	double offset_y = atof(pick_option(&c, &v, "-offset_y", "0"));
	double offset_z = atof(pick_option(&c, &v, "-offset_z", "0"));
	if (c != 4)
		return fprintf(stderr, "usage:\n\t"
			"%s dsm.tif zone rpc \n",*v);
			//0 1       2    4   5
	char *filename_dsm = v[1];
	int signed_zone    = atoi(v[2]);
	char *filename_rpc = v[3];

	// variables that will hold the local georeferencing transform
	double origin[2] = {0, 0};
	double scale[2] = {1, 1};

	// read georeferencing transform using GDAL
	GDALAllRegister();
	GDALDatasetH gdal_dataset = GDALOpen(filename_dsm, GA_ReadOnly);
	if (gdal_dataset == NULL)
		return fprintf(stderr, "GDALOpen(%s) failed\n", filename_dsm);
	double tmp[6];
	if (GDALGetGeoTransform(gdal_dataset, tmp) == CE_None) {
		origin[0] = tmp[0], origin[1] = tmp[3];
		scale[0] = tmp[1], scale[1] = tmp[5];
	} else {
		fprintf(stderr, "WARNING: not found origin and scale info\n");
	}

	// read the whole input DSM (typically, rather small)
	int w, h;
	float *x = iio_read_image_float(filename_dsm, &w, &h);
	if (!x)
		return fprintf(stderr, "iio_read(%s) failed\n", filename_dsm);

        // get mean_value of x
        float mean = 0;
        int count = 0;
        for (int i = 0; i < w*h; i++)
                if (x[i] > -9000)
                {
                        count +=1;
                        mean += x[i];
                }
        mean /= count;

	// open the reference rpc
	struct rpc huge_rpc[1];
	read_rpc_file_xml(huge_rpc, filename_rpc);

        // create P matrix and fill its last line
        double P[12];
        for (int i = 8; i < 11; i++)
                P[i]=0;
        P[11] = 1;

        // fill last column of P
        int i_centre = scale[0] * w/2 + origin[0];
        int j_centre = scale[1] * h/2 + origin[1];

        double ij[2];
        double z_centre = 20;

        float incrementation[12] = {0};
        for (int i = 0; i < 3; i++)
                incrementation[i*4] = 1;

        //incrementation[0] = scale[0];
        //incrementation[4] = scale[1];
        //incrementation[8] = 1;


        // evaluate projection for centre and 3 pixels in x, y, and z directions
        // and fill P with result
        for (int i = 0; i < 4; i++)
        {
                eval_proj(ij, i_centre + incrementation[i*3], 
                                j_centre + incrementation[i*3+1], 
                                z_centre + incrementation[i*3+2], 
                                huge_rpc, signed_zone);
                for (int j = 0; j < 2; j++)
                        P[j*4+i] = ij[j];
        }

        // get finite difference for first 3 columns
        for (int i = 0; i < 3; i++)
                for (int j = 0; j < 2; j++)
                        P[j*4+i] -= P[j*4+3];

       // for (int i = 0; i < 2; i++)
       //         for (int j = 0; j < 2; j++)
       //                 P[j*4+i] /= scale[i];

        // compute last column
        for (int j = 0; j < 2; j++)
                P[j*4+3] = P[j*4+3] - i_centre*P[j*4] - j_centre*P[j*4+1] 
                        - z_centre*P[j*4+2];


        // print P
        for (int j = 0; j < 3; j++)
        {
                for (int i = 0; i < 4; i++)
                        printf("%.20lf ",P[j*4+i]);
                printf("\n");
        }

        double z = 120;

        float *diff = malloc(2 * w * h * sizeof(float));
        for (int j=0; j<h; j++)
        for (int i=0; i<w; i++)
        {
                double ij_exact[2];
                double ij_approx[3];
                double e = i * scale[0] + origin[0];
                double n = j * scale[1] + origin[1];
                double xyz1[4] = {e, n, z, 1};

                eval_proj(ij_exact, e, n, z, 
                                huge_rpc, signed_zone);

                matrix_product_mxn_nxp(ij_approx, P, xyz1, 3, 4, 1);

                diff[j*w+i] = ij_exact[0]-ij_approx[0];
                diff[w*h+j*w+i] = ij_exact[1]-ij_approx[1]; 
        }
        iio_save_image_float_split("errors.tif", diff, w, h, 2);

        return 0;
}

int main(int c, char *v[]) { return main_get_projection_matrix(c,v); }
