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

void sun_plan_projection(double ijh[3], double eastnorthheight[3], 
        double scale[2], double az, double el)
{
    double e = eastnorthheight[0]; 
    double n = eastnorthheight[1]; 
    double z = eastnorthheight[2]; 

    ijh[0] = (e + (ijh[2] - z)*cos(el * M_PI /180) * sin(az * M_PI /180)
            /sin(el * M_PI /180))/scale[0];
    ijh[1] = -(n + (ijh[2] - z)*cos(el * M_PI /180) * cos(az * M_PI /180)
            /sin(el * M_PI /180))/scale[1];
}
// main_project: extracts corners position of lidar projection on huge image
//
// input_1 : a georeferenced DSM of size WxH with nan when no data
// input_2 : utm zone
// input_3 : the RPC of a georeferenced d-band image (of any size, typically huge) 
// input_4 : the RPC of the msi image
// output  : print xmin ymin W H rectangle to crop from the huge image
//
// note: formally, for the input we have to give the signed zone
//
int main_get_corners(int c, char *v[])
{
        double offset_x = atof(pick_option(&c, &v, "-offset_x", "0"));
        double offset_y = atof(pick_option(&c, &v, "-offset_y", "0"));
        double offset_z = atof(pick_option(&c, &v, "-offset_z", "0"));
        double bias_x = atof(pick_option(&c, &v, "bx", "0"));
        double bias_y = atof(pick_option(&c, &v, "by", "0"));
        if (c < 4)
                return fprintf(stderr, "usage:\n\t" 
                        "%s dsm_nan.tif zone rpc \n",*v);
                        //0 1           2    3
        char *filename_dsm = v[1];
        int signed_zone    = atoi(v[2]);
        double azimuth = atof(v[3]);
        double elevation = atof(v[4]);


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

        // get min and max height in lidar
        float zm[2] = {0, 0};

        for (int i = 0; i < w; i++)
                for (int j = 0; j < h; j++)
                {
                        if (x[i+j*w]>zm[1])
                                zm[1] = x[i+j*w];
                        if (x[i+j*w]<zm[0] && x[i+j*w] > -9000)
                                zm[0] = x[i+j*w];
                }
        zm[0] -= 2;
        zm[1] += 2;

        // initialize mins and maxs
        int xmin = 100000000;
        int ymin = 100000000;
        int xmax = -100000000;
        int ymax = -100000000;

        // for loop on corners of 3d dsm image to get projections on huge_image
        for (int l = 0; l < 2; l++)
        for (int j = 0; j < 2; j++)
        for (int i = 0; i < 2; i++)
        {
                if (!isfinite(x[j*(h-1)*w+i*(w-1)])) continue;
                double e = i*(w-1) * scale[0]; // easting
                double n = j*(h-1) * scale[1]; // northing
                double z = zm[l]; // height estimation

                // compute coordinates in huge image
                double eastnorthheight[3] = {e, n, z};
                double ijh[3] = {0, 0, 10000};
                sun_plan_projection(ijh, eastnorthheight, scale,
                        azimuth, elevation);

                // compare to mins and maxs
                if (ijh[0]<xmin)
                    xmin = (int) floor(ijh[0]);
                if (ijh[0]>xmax)
                    xmax = (int) ceil(ijh[0]);
                if (ijh[1]<ymin)
                    ymin = (int) floor(ijh[1]);
                if (ijh[1]>ymax)
                    ymax = (int) ceil(ijh[1]);

        }

                        
        printf("%d %d %d %d\n", xmin, ymin, xmax-xmin, ymax-ymin);

        // save and exit without cleanup

        return 0;
}

int main(int c, char *v[]) { return main_get_corners(c,v); }
