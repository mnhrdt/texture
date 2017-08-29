#include <assert.h>
#include <stdio.h>
#include <gdal.h>
#include <cpl_conv.h>
#include <libgen.h>

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

bool xy_are_in_bounds(float *xy, // tableau de coords (x1, x2,..., y1,y2,..)
                int l, // longueur des x/y
                float xmin, float xmax, float ymin, float ymax)
{       
        bool a = true;
        for (int i = 0; i < l; i++) 
                a &= (xy[i] >= xmin && xy[i] <= xmax &&
                                xy[i+l] >= ymin && xy[i+l] <=ymax);
        return a;
}


struct image_coord{ // coordonnées de la projection du sommet sur l'image
        float i; 
        float j;
};

struct vertex{
        float x; // x coord dans l'espace
        float y; // y coord dans l'espace
        float z; // z coord dans l'espace
        struct image_coord *im; // liste coordonnées sommet dans chaque im 
};

struct face{
        int v0; // 1er sommet
        int v1; // 2ème sommet
        int v2; // 3ème sommet
        float n[3]; // vecteur normal à la face
        int im; // référence de l'image la plus en face
};

struct mesh_t{
        int ni; // nombre d'images
        int nv; // nombre de sommets
        int nf; // nombre de faces
        struct vertex *v; // liste des sommets
        struct face *f; // liste des faces
};

bool ith_face_is_visible(struct mesh_t mesh, int i)
{
        struct face mf = mesh.f[i];
        return !isnan(mesh.v[mf.v0].im[mf.im].i) && 
                !isnan(mesh.v[mf.v0].im[mf.im].j) && 
                !isnan(mesh.v[mf.v1].im[mf.im].i) && 
                !isnan(mesh.v[mf.v1].im[mf.im].j) && 
                !isnan(mesh.v[mf.v2].im[mf.im].i) &&
                !isnan(mesh.v[mf.v2].im[mf.im].j);
}

int write_ply_t(char *filename_ply, char *filename_a, struct mesh_t mesh)
{
	// dump the ply file (with dsm-inherited connectivity)
	FILE *f = fopen(filename_ply, "w");
	if (!f) printf("WARNING: couldn't open ply file.\n");
        
        // print header
	fprintf(f, "ply\n");
	fprintf(f, "format ascii 1.0\n");
	fprintf(f, "comment created by cutrecombine\n");
        fprintf(f, "comment TextureFile %s.png\n", basename(filename_a)); // ici le plyfile note le nom de l'atlas avec le chemin d'accès. Je voudrais juste le nom de fichier, sans chemin d'accès avant.
	// if (offset_x) fprintf(f, "comment offset_x = %lf\n", offset_x);
	// if (offset_y) fprintf(f, "comment offset_y = %lf\n", offset_y);
	// if (offset_z) fprintf(f, "comment offset_z = %lf\n", offset_z);
	fprintf(f, "element vertex %d\n", mesh.nv);
	fprintf(f, "property float x\n");
	fprintf(f, "property float y\n");
	fprintf(f, "property float z\n");
	fprintf(f, "element face %d\n", mesh.nf);
	fprintf(f, "property list uchar int vertex_indices\n");
        fprintf(f, "property list uchar float texcoord\n");
	fprintf(f, "end_header\n");

        // print vertices
        for (int i = 0; i < mesh.nv; i++)
		fprintf(f, "%.16lf %.16lf %.16lf\n", 
                                mesh.v[i].x, mesh.v[i].y, mesh.v[i].z);

        // print faces
        for (int i = 0; i < mesh.nf; i++) {
                struct face mf = mesh.f[i];
                float a[6] = {mesh.v[mf.v0].im[mf.im].i, mesh.v[mf.v1].im[mf.im].i,
                        mesh.v[mf.v2].im[mf.im].i, mesh.v[mf.v0].im[mf.im].j,
                        mesh.v[mf.v1].im[mf.im].j, mesh.v[mf.v2].im[mf.im].j};
                if (ith_face_is_visible(mesh, i))
                {
                        fprintf(f, "3 %d %d %d 6 %f %f %f %f %f %f \n", 
                                mf.v0, mf.v1, mf.v2, 
                                a[0], -a[3]/2, a[1], -a[4]/2, a[2], -a[5]/2);
                        if (!xy_are_in_bounds(a, 3, 0, 1, 0, 1))
                                return fprintf(stderr, "WARNING: tries to access invalid texture %f %f %f %f %f %f\n",a[0], a[3], a[1], a[4], a[2], a[5]);
                }
                else
                        fprintf(f, "3 %d %d %d 6 0 0 0 0.1 0.1 0\n",
                                mf.v0, mf.v1, mf.v2);

        }
        return 0;
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
        double ox = atof(pick_option(&c, &v, "ox", "0"));
        double oy = atof(pick_option(&c, &v, "oy", "0"));
	if (c < 6)
		return fprintf(stderr, "usage:\n\t"
			"%s dsm.tif img.tif match.tif .ply atlas\n",*v);
			//0 1       2       3         4    5
	char *filename_dsm = v[1];
	char *filename_img = v[2];
	char *filename_m = v[3];
	char *filename_ply = v[4];
        char *filename_a = v[5];
        printf("ox %f oy %f\n", ox, oy);

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

        printf("scale x %f scale y %f\n", scale[0], scale[1]);
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

        // create mesh
        struct mesh_t mesh;
        mesh.ni = 1;

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
        mesh.nv = nvertices;

	// count number of valid square faces
	int nfaces = 0;
	for (int j = 0; j < h-1; j++)
	for (int i = 0; i < w-1; i++)
	{
		int q[4] = {vid[j][i], vid[j+1][i], vid[j+1][i+1], vid[j][i+1]};
		if (q[0] >= 0 && q[1] >= 0 && q[2] >= 0 && q[3] >= 0)
			nfaces += 2;
	}
        mesh.nf = nfaces;

	// output vertices
        mesh.v = malloc(mesh.nv*sizeof(struct vertex));
	int cx = 0;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		if (!isfinite(height[j][i])) continue;
		// compute lonlat from eastnorth = {p[0], p[1]}
		double e = i * scale[0];// + origin[0]; // easting
		double n = j * scale[1];// + origin[1]; // northing
		double z = height[j][i];             // height
                mesh.v[cx].x = e - offset_x;
                mesh.v[cx].y = n - offset_y;
                mesh.v[cx].z = z - offset_z;

                mesh.v[cx].im = malloc((mesh.ni)*sizeof(struct image_coord));
                mesh.v[cx].im[0].i = (m[2*(i+j*w)]+ox)  /wi;
                mesh.v[cx].im[0].j = (m[2*(i+j*w)+1]+oy)/hi;
		cx += 1;
	}
        printf("offset x %f y %f\n", ox, oy);
	assert(cx == nvertices);

	// output faces
        mesh.f = malloc(mesh.nf*sizeof(struct face));
	cx = 0;
	for (int j = 0; j < h-1; j++)
	for (int i = 0; i < w-1; i++)
	{
		int q[4] = {vid[j][i], vid[j+1][i], vid[j+1][i+1], vid[j][i+1]};
		if (q[0] >= 0 && q[1] >= 0 && q[2] >= 0 && q[3] >= 0)
		{
                        mesh.f[cx].im = 0;
                        mesh.f[cx] = (struct face) {.v0 = q[3], 
                                .v1 = q[0], .v2 = q[1]};
                        mesh.f[cx+1] = (struct face) {.v0 = q[3], 
                                .v1 = q[1], .v2 = q[2]};
			cx += 2;
		}
	}
	assert(cx == nfaces);

        write_ply_t(filename_ply, filename_a, mesh);
	return 0;
}

int main(int c, char *v[]) { return main_colorsingle(c,v); }
