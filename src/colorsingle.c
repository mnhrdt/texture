#include <assert.h>
#include <stdio.h>
#include <gdal.h>
#include <cpl_conv.h>

#include "rpc.c"
#include "fail.c"
#include "iio.h"

// extrapolate by nearest value
static float getpixel_f(float *I, int w, int h, int i, int j)
{
        if (i < 0) i = 0;
        if (j < 0) j = 0;
        if (i >= w) i = w-1;
        if (j >= h) j = h-1;
        return I[i+j*w];
}

// elevate a georeferenced DSM to an ascii point cloud
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
int main_colorsingle(int c, char *v[])
{
	double offset_x = atof(pick_option(&c, &v, "-offset_x", "0"));
	double offset_y = atof(pick_option(&c, &v, "-offset_y", "0"));
	double offset_z = atof(pick_option(&c, &v, "-offset_z", "0"));
	if (c != 6)
		return fprintf(stderr, "usage:\n\t"
			"%s dsm.tif img.tif match.tif .ply atlas\n",*v);
			//0 1       2       3         4    5
	char *filename_dsm = v[1];
	char *filename_img = v[2];
	char *filename_m = v[3];
	char *filename_ply = v[4];
        char *filename_a = v[5];

        // variables that will hold the local georeferencing transform
        double origin[2] = {0, 0};
        double scale[2] = {1, 1};

	// read the whole input DSM (typically, rather small)
	int w, h;
	float *x = iio_read_image_float(filename_dsm, &w, &h);
	if (!x)
		return fprintf(stderr, "iio_read(%s) failed\n", filename_dsm);

        // read lidar projections from colorize_with_shadows
        int wm, hm, pdm;
        float *m = iio_read_image_float_vec(filename_m, &wm, &hm, &pdm);
        if (!m)
                return fprintf(stderr, "iio_read(%s) failed\n", filename_m);
        if (w != wm || h != hm || pdm != 2)
                return fprintf(stderr, "input sizes mismatch (%s)\n", filename_m);

        // read the image to use as texture (crop corresponding to the dsm)
        int wi, hi;
        float *img = iio_read_image_float(filename_img, &wi, &hi);
        if (!img)
                return fprintf(stderr, "iio_read(%s) failed\n", filename_img);

        // create and save atlas by adding bright green at the bottom at the image
        float *a = malloc(6 * wi * hi * sizeof(float));
        for (int j = 0; j < hi; j++)
        for (int i = 0; i < wi; i++)
        for (int l = 0; l < 3; l++)
        {
                a[(i+j*wi)*3+l] = getpixel_f(img, wi, hi, i, j); 
                a[(i+j*wi)*3+l+wi*hi*3] = (l==1) ? 2055 : 0;
        }

        char n_a[1000];
        sprintf(n_a,"%s.tif",filename_a);
        iio_save_image_float_vec(n_a, a, wi, 2*hi, 3);


	// dump the ply file (with dsm-inherited connectivity)
	FILE *f = fopen(filename_ply, "w");
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

	// count number of valid square faces
	int nfaces = 0;
	for (int j = 0; j < h-1; j++)
	for (int i = 0; i < w-1; i++)
	{
		int q[4] = {vid[j][i], vid[j+1][i], vid[j+1][i+1], vid[j][i+1]};
		if (q[0] >= 0 && q[1] >= 0 && q[2] >= 0 && q[3] >= 0)
			nfaces += 2;
	}

	// print header
	fprintf(f, "ply\n");
	fprintf(f, "format ascii 1.0\n");
	fprintf(f, "comment created by cutrecombine\n");
        fprintf(f, "comment TextureFile %s.png\n", filename_a); // ici le plyfile note le nom de l'atlas avec le chemin d'accès. Je voudrais juste le nom de fichier, sans chemin d'accès avant.
	if (offset_x) fprintf(f, "comment offset_x = %lf\n", offset_x);
	if (offset_y) fprintf(f, "comment offset_y = %lf\n", offset_y);
	if (offset_z) fprintf(f, "comment offset_z = %lf\n", offset_z);
	fprintf(f, "element vertex %d\n", nvertices);
	fprintf(f, "property float x\n");
	fprintf(f, "property float y\n");
	fprintf(f, "property float z\n");
	fprintf(f, "element face %d\n", nfaces);
	fprintf(f, "property list uchar int vertex_indices\n");
        fprintf(f, "property list uchar float texcoord\n");
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
		double xyz[3] = {e - offset_x, n - offset_y, z - offset_z };
		fprintf(f, "%.16lf %.16lf %.16lf\n",
				xyz[0], xyz[1], xyz[2]);
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
                        if (m[2*(i+j*w)] == NAN || m[2*(i+(j+1)*w)] == NAN || 
                                        m[2*(i+1+(j+1)*w)] == NAN || 
                                        m[2*(i+1+j*w)] == NAN)
                        {
                                fprintf(f, "3 %d %d %d 6 %f %f %f %f %f %f \n", 
                                        q[3], q[1], q[0], 
                                        m[2*(i+1+ j   *w)]/wi, m[2*(i+1+ j   *w)+1]/(2*hi),
                                        m[2*(i+  (j+1)*w)]/wi, m[2*(i+  (j+1)*w)+1]/(2*hi),
                                        m[2*(i+   j   *w)]/wi, m[2*(i+   j   *w)+1]/(2*hi));
			        fprintf(f, "3 %d %d %d 6 %f %f %f %f %f %f \n", 
                                        q[3], q[2], q[1], 
                                        m[2*(i+1+ j   *w)]/wi, m[2*(i+1+ j   *w)+1]/(2*hi),
                                        m[2*(i+1+(j+1)*w)]/wi, m[2*(i+1+(j+1)*w)+1]/(2*hi),
                                        m[2*(i+  (j+1)*w)]/wi, m[2*(i+  (j+1)*w)+1]/(2*hi));
                        }
                        else
                        {
                                fprintf(f, "3 %d %d %d 6 %f %f %f %f %f %f \n", 
                                        q[3], q[1], q[0], 
                                        m[2*(i+1+ j   *w)]/wi, -m[2*(i+1+ j   *w)+1]/(2*hi),
                                        m[2*(i+  (j+1)*w)]/wi, -m[2*(i+  (j+1)*w)+1]/(2*hi),
                                        m[2*(i+   j   *w)]/wi, -m[2*(i+   j   *w)+1]/(2*hi));
			        fprintf(f, "3 %d %d %d 6 %f %f %f %f %f %f \n", 
                                        q[3], q[2], q[1], 
                                        m[2*(i+1+ j   *w)]/wi, -m[2*(i+1+ j   *w)+1]/(2*hi),
                                        m[2*(i+1+(j+1)*w)]/wi, -m[2*(i+1+(j+1)*w)+1]/(2*hi),
                                        m[2*(i+  (j+1)*w)]/wi, -m[2*(i+  (j+1)*w)+1]/(2*hi));
                        }
			cx += 2;
		}
	}
	assert(cx == nfaces);
	return 0;
}

int main(int c, char *v[]) { return main_colorsingle(c,v); }
