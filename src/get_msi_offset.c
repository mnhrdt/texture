#include <assert.h>
#include <stdio.h>
#include <gdal.h>
#include <cpl_conv.h>

#include "rpc.c"
#include "fail.c"
#include "iio.h"
#include "pickopt.c"

int utm_from_lonlat(double out_eastnorth[2], double lon, double lat);
void lonlat_from_eastnorthzone(double out_lonlat[2], double e, double n, int z);


// main_project: extracts corners position of lidar projection on huge image
//
// input_1 : a georeferenced DSM of size WxH with nan when no data
// input_2 : utm zone
// input_3 : the RPC of a georeferenced d-band image (of any size, typically huge) 
// output  : print xmin ymin W H rectangle to crop from the huge image
//
// note: formally, for the input we have to give the signed zone
//
int main_get_msi_offset(int c, char *v[])
{
        if (c < 5)
                return fprintf(stderr, "usage:\n\t"
                                "%s dsm_nan.tif zone pan.rpc msi.xml\n",*v);
                                //0 1           2    3       4
        char *filename_dsm = v[1];
        int signed_zone    = atoi(v[2]);
        char *filename_rpc = v[3];
        char *filename_xml = v[4];


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


        // open the reference rpc
        struct rpc huge_rpc[1];
        read_rpc_file_xml(huge_rpc, filename_rpc);

        struct rpc huge_xml[1];
        read_rpc_file_xml(huge_xml, filename_xml);
        double e = w/2 * scale[0] + origin[0];
        double n = h/2 * scale[1] + origin[1];
        double z = x[(w/2)+w*h/2]+20;

        double lonlat[3] = {0, 0, z};
        lonlat_from_eastnorthzone(lonlat, e, n, signed_zone);

        // compute coordinates in huge image
        double ij[2];
        rpc_projection(ij, huge_rpc, lonlat);
        double ij2[2];
        rpc_projection(ij2, huge_xml, lonlat);
        printf("%lf %lf\n", ij2[0]-ij[0]/4, ij2[1]-ij[1]/4);

                        

        // save and exit without cleanup

        return 0;
}

int main(int c, char *v[]) { return main_get_msi_offset(c,v); }
