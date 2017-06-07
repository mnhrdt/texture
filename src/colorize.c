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
	int r = GDALRasterIO(img, GF_Read, i,j,1, 1, roi,1,1, GDT_Float32, 0,0);
	return roi[0*0+0];
}

// elevate a georeferenced DSM to an ascii point cloud
// optionally, colorize the point cloud from a given reference image
// optionally, create a ply file with the desired connectivity
#include "pickopt.c"
int main_elevate(int c, char *v[])
{
#if 0
	char *out_ply = pick_option(&c, &v, "p", "");
	double offset_x = atof(pick_option(&c, &v, "-offset_x", "0"));
	double offset_y = atof(pick_option(&c, &v, "-offset_y", "0"));
	double offset_z = atof(pick_option(&c, &v, "-offset_z", "0"));
	if (c != 5 && c != 2)
		return fprintf(stderr, "usage:\n\t"
			"%s dsm [zone img rpc] [-p out.ply] > xyzi.txt\n",*v);
			//0 1    2    3   4
	char *filename_dsm = v[1];
	int signed_zone    = c > 2 ? atoi(v[2]) : 0;
	char *filename_img = c > 2 ? v[3] : NULL;
	char *filename_rpc = c > 2 ? v[4] : NULL;

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

	// if necessary, open the reference image and its rpc
	GDALRasterBandH huge_img[3];
	struct rpc huge_rpc[1];
	int pd = 1;
	if (filename_img) {
		read_rpc_file_xml(huge_rpc, filename_rpc);
		GDALDatasetH huge_dataset = GDALOpen(filename_img, GA_ReadOnly);
		pd = GDALGetRasterCount(huge_dataset);
		for (int i = 0; i < pd; i++)
			huge_img[i] = GDALGetRasterBand(huge_dataset, i+1);
	}

	// elevate each non-nan pixel
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		if (!isfinite(x[j*w+i])) continue;
		double e = i * scale[0] + origin[0]; // easting
		double n = j * scale[1] + origin[1]; // northing
		double z = x[j*w+i];                 // height
		float intensity = 0;
		if (filename_img) {
			// compute lonlat from eastnorth = {p[0], p[1]}
			double lonlat[3] = {0, 0, z};
			lonlat_from_eastnorthzone(lonlat, e, n, signed_zone);

			// compute coordinates in huge image
			double ij[2];
			rpc_projection(ij, huge_rpc, lonlat);

			// evaluate the color at this point
			intensity = gdal_getpixel(huge_img[0], ij[0], ij[1]);
		}
		e = e - offset_x;
		n = n - offset_y;
		z = z - offset_z;
		printf("%a %a %lf %g\n", e, n, z, intensity);
	}

	// if no ply file is requested, exit now
	if ('\0' == *out_ply)
		return 0;

	// dump the ply file (with dsm-inherited connectivity)
	FILE *f = fopen(out_ply, "w");
	if (!f) return 1;

	// assign comfortable pointers
	float (*height)[w] = (void*)x;
	int (*vid)[w] = malloc(w*h*sizeof(int));

	// count number of valid vertices
	int nvertices = 0;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
		if (isfinite(height[j][i]))
			vid[j][i] = nvertices++;
		else
			vid[j][i] = -1;

	// count number of valid faces
	int nfaces = 0;
	for (int j = 0; j < h-1; j++)
	for (int i = 0; i < w-1; i++)
	{
		int q[4] = {vid[j][i], vid[j+1][i], vid[j+1][i+1], vid[j][i+1]};
		if (q[0] >= 0 && q[1] >= 0 && q[2] >= 0 && q[3] >= 0)
			nfaces += 1;
	}

	// print header
	fprintf(f, "ply\n");
	fprintf(f, "format ascii 1.0\n");
	fprintf(f, "comment created by cutrecombine\n");
	if (offset_x) fprintf(f, "comment offset_x = %lf\n", offset_x);
	if (offset_y) fprintf(f, "comment offset_y = %lf\n", offset_y);
	if (offset_z) fprintf(f, "comment offset_z = %lf\n", offset_z);
	fprintf(f, "element vertex %d\n", nvertices);
	fprintf(f, "property float x\n");
	fprintf(f, "property float y\n");
	fprintf(f, "property float z\n");
	fprintf(f, "property uchar red\n");
	fprintf(f, "property uchar green\n");
	fprintf(f, "property uchar blue\n");
	fprintf(f, "element face %d\n", nfaces);
	fprintf(f, "property list uchar int vertex_index\n");
	fprintf(f, "end_header\n");
	int cx;

	// output vertices
	cx = 0;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		if (!isfinite(height[j][i])) continue;
		// compute lonlat from eastnorth = {p[0], p[1]}
		double e = i * scale[0] + origin[0]; // easting
		double n = j * scale[1] + origin[1]; // northing
		double z = height[j][i];             // height
		float color[3] = {200, 200, 200};
		if (filename_img) {
			double lonlat[3] = {0, 0, z};
			lonlat_from_eastnorthzone(lonlat, e, n, signed_zone);

			// compute coordinates in huge image
			double ij[2];
			rpc_projection(ij, huge_rpc, lonlat);

			// evaluate the color at this point
			for (int k = 0; k < pd; k++)
				color[k] = gdal_getpixel(huge_img[k], ij[0], ij[1]);
			for (int k = pd; k < 3; k++) color[k] = color[k-1];
		}
		uint8_t rgb[3] = { color[0], color[1], color[2] };
		double xyz[3] = {e - offset_x, n - offset_y, z - offset_z };
		fprintf(f, "%.16lf %.16lf %.16lf %d %d %d\n",
				xyz[0], xyz[1], xyz[2], rgb[0], rgb[1], rgb[2]);
		cx += 1;
	}
	assert(cx == nvertices);

	// output faces
	cx = 0;
	for (int j = 0; j < h-1; j++)
	for (int i = 0; i < w-1; i++)
	{
		int q[4] = {vid[j][i], vid[j+1][i], vid[j+1][i+1], vid[j][i+1]};
		if (q[0] >= 0 && q[1] >= 0 && q[2] >= 0 && q[3] >= 0)
		{
			fprintf(f, "4 %d %d %d %d\n", q[0], q[1], q[2], q[3]);
			cx += 1;
		}
	}
	assert(cx == nfaces);
#endif
	return 0;
}

// main_colorize: like elevate, but produce a color image
//
// input_1 : a georeferenced DSM of size WxH
// input_2 : a georeferenced d-band image (of any size, typically huge)
// input_3 : the RPC of the d-band image
// output  : a d-band image of size WxH, with colors from the correct place
//
// note: formally, for the input we have to give the signed zone
//
int main_colorize(int c, char *v[])
{
	double offset_x = atof(pick_option(&c, &v, "-offset_x", "0"));
	double offset_y = atof(pick_option(&c, &v, "-offset_y", "0"));
	double offset_z = atof(pick_option(&c, &v, "-offset_z", "0"));
	if (c != 6)
		return fprintf(stderr, "usage:\n\t"
			"%s dsm.tif zone img.ntf rpc out.tif\n",*v);
			//0 1       2    3       4   5
	char *filename_dsm = v[1];
	int signed_zone    = atoi(v[2]);
	char *filename_img = v[3];
	char *filename_rpc = v[4];
	char *filename_out = v[5];

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

	// open the reference image and obtain its pixel dimension "pd"
	GDALDatasetH huge_dataset = GDALOpen(filename_img, GA_ReadOnly);
	int pd = GDALGetRasterCount(huge_dataset);
	GDALRasterBandH huge_img[pd];
	for (int i = 0; i < pd; i++)
		huge_img[i] = GDALGetRasterBand(huge_dataset, i+1);

	// open the reference rpc
	struct rpc huge_rpc[1];
	read_rpc_file_xml(huge_rpc, filename_rpc);

	// allocate space for output image
	float *out_x = malloc(w * h * pd * sizeof*out_x);

	// elevate each non-nan pixel
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		if (!isfinite(x[j*w+i])) continue;
		double e = i * scale[0] + origin[0]; // easting
		double n = j * scale[1] + origin[1]; // northing
		double z = x[j*w+i];                 // height
		float intensity = 0;
		float pix[pd];
		if (filename_img) {
			// compute lonlat from eastnorth = {p[0], p[1]}
			double lonlat[3] = {0, 0, z};
			lonlat_from_eastnorthzone(lonlat, e, n, signed_zone);

			// compute coordinates in huge image
			double ij[2];
			rpc_projection(ij, huge_rpc, lonlat);

			// evaluate the color at this point
			for (int l = 0; l < pd; l++)
				pix[l] = gdal_getpixel(huge_img[l],ij[0],ij[1]);
		}
		// TODO: decide whether the offset is added to the image or
		// to the DSM (probably the image, but I'm not sure yet)
		//e = e - offset_x;
		//n = n - offset_y;
		//z = z - offset_z;
		for (int l = 0; l < pd; l++)
			out_x[(j*w+i)*pd + l] = pix[l];
	}

	// save and exit without cleanup
	iio_save_image_float_vec(filename_out, out_x, w, h, pd);
	return 0;
}

int main(int c, char *v[]) { return main_colorize(c,v); }
