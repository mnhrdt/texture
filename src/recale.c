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

int utm_from_lonlat(double out_eastnorth[2], double lon, double lat);
void lonlat_from_eastnorthzone(double out_lonlat[2], double e, double n, int z);




#include "pickopt.c"

int main_recale(int c, char *v[])
{
        double offset_x = atof(pick_option(&c, &v, "-offset_x", "0"));
	double offset_y = atof(pick_option(&c, &v, "-offset_y", "0"));
	double offset_z = atof(pick_option(&c, &v, "-offset_z", "0"));
	if (c < 4)
		return fprintf(stderr, "usage:\n\t"
			"%s dsm.tif rpc_i.tif ncc.txt\n",*v);
			//0 1       2         3      
	char *filename_dsm = v[1];
        char *filename_rpc = v[2];
	char *filename_ncc = v[3];
        
        // variables that will hold the local georeferencing transform
        double origin[2] = {0, 0};
        double scale[2] = {1, 1};

        // read georeferencing transform using GDAL
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

        FILE *bias;
        bias = fopen(filename_ncc,"r");
        if (!bias)
                return fprintf(stderr, "fopen(%s) failed\n", filename_ncc);
        double xysz[4];
        for (int i = 0; i < 4; i++)
        {   
                int r = fscanf(bias, "%lf", &xysz[i]);
                if (r!=1)
                        return fprintf(stderr, "could not read element %d of %s", i, filename_ncc);
        }  
        double offset[3] = {xysz[0]*scale[0], xysz[1]*scale[1], xysz[3]};

        struct rpc huge_rpc[1];
        read_rpc_file_xml(huge_rpc, filename_rpc);
        

        double z = huge_rpc[0].offset[2];
        double ijh[3] = {huge_rpc[0].offset[0],huge_rpc[0].offset[1],huge_rpc[0].offset[2]};
        printf("i %lf j %lf z %lf\n", ijh[0], ijh[1], z);
        double lonlat[2] = {0, 0};
        rpc_localization(lonlat, huge_rpc, ijh);
        double en[2];
        utm_from_lonlat(en, lonlat[0], lonlat[1]);
        printf("e %lf n %lf\n", en[0], en[1]);
        for (int i = 0; i < 2; i++)
                en[i] += offset[i];
        printf("e %lf n %lf z %lf\n", en[0], en[1], z+offset[2]);
        lonlat_from_eastnorthzone(lonlat, en[0], en[1], -21);

        double ij[2] = {0, 0};
        double lonlatheight[3] = {lonlat[0], lonlat[1], z+offset[2]};
        rpc_projection(ij, huge_rpc, lonlatheight);
        printf("i %lf j %lf z %lf\n", ij[0], ij[1], z);

        printf("décalage théorique en pixel %lf %lf\n", ij[0]-ijh[0], ij[1]-ijh[1]);


        double ijh2[3] = {ijh[0]-4, ijh[1]+10, z};
        rpc_localization(lonlat, huge_rpc, ijh2);
        double en2[2]; 
        utm_from_lonlat(en2, lonlat[0], lonlat[1]);
	return 0;
        
}

int main(int c, char *v[]) { return main_recale(c,v); }

