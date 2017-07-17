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

// extrapolate by nearest value
static float getpixel_f(float *I, int w, int h, int i, int j)
{
        if (i < 0) i = 0;
        if (j < 0) j = 0;
        if (i >= w) i = w-1;
        if (j >= h) j = h-1;
        return I[i+j*w];
}

double scalar_product(double *u, double *v, int dim)
{
        double sp = 0;
        for (int i = 0; i < dim; i++)
                sp += u[i]*v[i];
        return sp;
}

void cross_product(double axb[3], double a[3], double b[3])
{
	// a0 a1 a2
	// b0 b1 b2
	axb[0] = a[1]*b[2] - a[2]*b[1];
	axb[1] = a[2]*b[0] - a[0]*b[2];
	axb[2] = a[0]*b[1] - a[1]*b[0];
}

double euclidean_norm(double *x, int n)
{
	return n > 0 ? hypot(*x, euclidean_norm(x+1, n-1)) : 0;
}


void triangle_normal(double n[3], double a[3], double b[3], double c[3]) // les sommets sont donnés dans le sens direct
{
        double u[3];
        double v[3];
        for (int i = 0; i < 3; i++) u[i] = b[i] - a[i];
        for (int i = 0; i < 3; i++) v[i] = c[i] - a[i];
	cross_product(n, u, v);

        double norm = euclidean_norm(n, 3);
        for (int i = 0; i < 3; i++) n[i] /= norm;

        norm = euclidean_norm(n, 3);
        if (norm < 0.9999 || norm > 1.0000001)
                printf("WARNING: normalisation error in triangle_normal, norme = %.16lf\n", norm);
}

bool xy_are_in_bounds(double *xy, // tableau de coords (x1, x2,..., y1,y2,..)
                int l, // longueur des x/y
                double xmin, double xmax, double ymin, double ymax)
{       
        bool a = true;
        for (int i = 0; i < l; i++) 
                a &= (xy[i] >= xmin && xy[i] <= xmax &&
                                xy[i+l] >= ymin && xy[i+l] <=ymax);
        return a;
}


struct image_coord{ // coordonnées de la projection du sommet sur l'image
        double i; 
        double j;
};

struct vertex{
        int ij[2]; // coordinates in lidar
        double xyz[3]; // coord dans l'espace
        struct image_coord *im; // liste coordonnées sommet dans chaque im 
};

struct face{
        int v0; // 1er sommet
        int v1; // 2ème sommet
        int v2; // 3ème sommet
        double n[3]; // vecteur normal à la face
        int im; // référence de l'image la plus en face
        bool *isvisible; // 1 si visible dans l'image i, 0 sinon
};

struct mesh_t{
        int nimages; // nombre d'images
        int nv; // nombre de sommets
        int nf; // nombre de faces
        struct vertex *v; // liste des sommets
        struct face *f; // liste des faces
};

bool ith_face_is_visible_in_image(struct mesh_t mesh, int i, int j)
{
        struct face mf = mesh.f[i];
        return !isnan(mesh.v[mf.v0].im[j].i) && 
                !isnan(mesh.v[mf.v0].im[j].j) && 
                !isnan(mesh.v[mf.v1].im[j].i) && 
                !isnan(mesh.v[mf.v1].im[j].j) && 
                !isnan(mesh.v[mf.v2].im[j].i) &&
                !isnan(mesh.v[mf.v2].im[j].j);
}

void camera_direction(double n[3], struct rpc *r)
{
        // initialise height and fill 3rd vector coordinate
        double z = 0;
        n[2] = -1;

        // get first 3D point using localisation 
        double ijh[3] = {500, 500, z};
        double lonlat[2] = {0, 0};
        rpc_localization(lonlat, r, ijh);
        double en[2];
        utm_from_lonlat(en, lonlat[0], lonlat[1]);
        for (int i = 0; i < 2; i++)
                n[i] = en[i];

        // get second 3D point combining localisation and projection
        double ij[2] = {0, 0};
        double lonlatheight[3] = {lonlat[0], lonlat[1], z};
        rpc_projection(ij, r, lonlatheight);
        for (int i = 0; i < 2; i++)
                ijh[i] = ij[i];
        ijh[2] = z + 1;
        rpc_localization(lonlat, r, ijh);
        utm_from_lonlat(en, lonlat[0], lonlat[1]);

        // fill in the first two vector coordinates
        for (int i = 0; i < 2; i++)
                n[i] -= en[i];

        // normalise direction vector
        for (int i = 0; i < 3; i++)
                n[i] /= sqrt((pow(n[0],2) + pow(n[1],2) + pow(n[2],2)));
}

void write_ply_t(char *filename_ply, char *filename_a, struct mesh_t mesh)
{
	// dump the ply file (with dsm-inherited connectivity)
	FILE *f = fopen(filename_ply, "w");
	if (!f) printf("WARNING: couldn't open ply file.\n");
        
        // print header
	fprintf(f, "ply\n");
	fprintf(f, "format ascii 1.0\n");
	fprintf(f, "comment created by cutrecombine\n");
        fprintf(f, "comment TextureFile %s.png\n", basename(filename_a));
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
                                mesh.v[i].xyz[0], mesh.v[i].xyz[1], mesh.v[i].xyz[2]);

        // print faces
        for (int i = 0; i < mesh.nf; i++) {
                struct face mf = mesh.f[i];
                double a[6] = {mesh.v[mf.v0].im[mf.im].i, mesh.v[mf.v1].im[mf.im].i,
                        mesh.v[mf.v2].im[mf.im].i, mesh.v[mf.v0].im[mf.im].j,
                        mesh.v[mf.v1].im[mf.im].j, mesh.v[mf.v2].im[mf.im].j};
                if (ith_face_is_visible_in_image(mesh, i, mf.im))
                {
                        fprintf(f, "3 %d %d %d 6 %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf \n", 
                                mf.v0, mf.v1, mf.v2, 
                                a[0], -a[3]/2, a[1], -a[4]/2, a[2], -a[5]/2);
                        if (!xy_are_in_bounds(a, 3, 0, 1, 0, 1))
                                printf("WARNING: tries to access invalid texture\n");
                }
                else
                        fprintf(f, "3 %d %d %d 6 0 0 0 0.1 0.1 0\n",
                                mf.v0, mf.v1, mf.v2);

        }
}

void write_ply_map_t(char *filename_ply, char *filename_a, struct mesh_t mesh)
{
	// dump the ply file (with dsm-inherited connectivity)
        char filename_map[1000];
        sprintf(filename_map, "%s/map.ply", dirname(filename_ply));
	FILE *f = fopen(filename_map, "w");
	if (!f) printf("WARNING: couldn't open ply file.\n");
        
        // print header
	fprintf(f, "ply\n");
	fprintf(f, "format ascii 1.0\n");
	fprintf(f, "comment created by cutrecombine\n");
        fprintf(f, "comment TextureFile %s_map.png\n", basename(filename_a));
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
                                mesh.v[i].xyz[0], mesh.v[i].xyz[1], mesh.v[i].xyz[2]);

        // print faces
        for (int i = 0; i < mesh.nf; i++) {
                struct face mf = mesh.f[i];
                double a[6] = {mesh.v[mf.v0].im[mf.im].i, mesh.v[mf.v1].im[mf.im].i,
                        mesh.v[mf.v2].im[mf.im].i, mesh.v[mf.v0].im[mf.im].j,
                        mesh.v[mf.v1].im[mf.im].j, mesh.v[mf.v2].im[mf.im].j};
                if (ith_face_is_visible_in_image(mesh, i, mf.im))
                {
                        fprintf(f, "3 %d %d %d 6 %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf \n", 
                                mf.v0, mf.v1, mf.v2, 
                                a[0], -a[3]/2, a[1], -a[4]/2, a[2], -a[5]/2);
                        if (!xy_are_in_bounds(a, 3, 0, 1, 0, 1))
                                printf("WARNING: tries to access invalid texture\n");
                }
                else
                        fprintf(f, "3 %d %d %d 6 0 0 0 0.1 0.1 0\n",
                                mf.v0, mf.v1, mf.v2);

        }
}

// elevate a georeferenced DSM to an ascii point cloud

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
int main_colormultiple(int c, char *v[])
{
        double offset_x = atof(pick_option(&c, &v, "-offset_x", "0"));
	double offset_y = atof(pick_option(&c, &v, "-offset_y", "0"));
	double offset_z = atof(pick_option(&c, &v, "-offset_z", "0"));
	if (c % 3 != 1 || c < 5)
		return fprintf(stderr, "usage:\n\t"
			"%s dsm.tif img_i.tif rpc_i.tif match_i.tif out.ply atlas\n",*v);
			//0 1       2         3         4           3*i+2   3*i+3
        int nimages = (c-4)/3;
	char *filename_dsm = v[1];
        char *filename_m = v[4];
	char *filename_ply = v[3*nimages+2];
        char *filename_a = v[3*nimages+3];

        // create and save atlas by putting all the images on a line, 
        // padding by green and adding green on the bottom half of the atlas.
        int wi, hi;
        int wimax = 0;
        int himax = 0;
        for (int ni = 0; ni < nimages; ni++)
        {
                char *filename_img = v[ni+2];
                float *img = iio_read_image_float(filename_img, &wi, &hi);
                if (!img)
                        return fprintf(stderr, "iio_read(%s) failed\n", filename_img);
                if (wi > wimax)
                        wimax = wi;
                if (hi > himax)
                        himax = hi;
                free(img);
        }
        int pda = 3;
        float *a = malloc(2 * pda * wimax * himax * nimages * sizeof(float));
        for (int i = 0; i < 2*wimax*himax*nimages; i++)
        for (int l = 0; l < pda; l++)
                a[3*i+l] = (l==1) ? 2055 : 0;
        for (int ni = 0; ni < nimages; ni++)
        {
                char *filename_img = v[ni+2];
                float *img = iio_read_image_float(filename_img, &wi, &hi);
                if (!img)
                        return fprintf(stderr, "iio_read(%s) failed\n", filename_img);
                for (int i = 0; i < wi; i++)
                for (int j = 0; j < hi; j++)
                for (int l = 0; l < pda; l++)
                        a[3*(i+ni*wimax+j*nimages*wimax)+l] = img[i+j*wi];
                free(img);
        }
        char n_a[1000];
        sprintf(n_a, "%s.tif", filename_a);
        iio_save_image_float_vec(n_a, a, nimages*wimax, 2*himax, pda);
        free(a);

        // get camera directions
        double *cam_n = malloc(3 * nimages * sizeof(double));
        double c_n[3] = {0};
        struct rpc huge_rpc[1];
        for (int ni = 0; ni < nimages; ni++)
        {
                char *filename_rpc = v[nimages+ni+2];
                read_rpc_file_xml(huge_rpc, filename_rpc);
                camera_direction(c_n, huge_rpc);
                for (int l = 0; l < 3; l++)
                        cam_n[3*ni+l] = c_n[l];
        }




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


        // create mesh
        struct mesh_t mesh;
        mesh.nimages = nimages;

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
		double e = i; //* scale[0] + origin[0]; // easting
		double n = j; //* scale[1] + origin[1]; // northing
		double z = height[j][i];             // height
                mesh.v[cx].xyz[0] = e - offset_x;
                mesh.v[cx].xyz[1] = n - offset_y;
                mesh.v[cx].xyz[2] = z - offset_z;

                mesh.v[cx].ij[0] = i;
                mesh.v[cx].ij[1] = j;

                mesh.v[cx].im = malloc((mesh.nimages)*sizeof(struct image_coord));
		cx += 1;
	}
	assert(cx == nvertices);

        // for each vertex, get texture coordinates
       int wm, hm, pdm;
       for (int ni = 0; ni < nimages; ni++)
       {
               char *filename_m = v[2*nimages+ni+2];
               float *m = iio_read_image_float_vec(filename_m, &wm, &hm, &pdm);
               if (!m)
                       return fprintf(stderr, "iio_read(%s) failed\n", filename_m);
               if (w != wm || h != hm || pdm != 2)
                       return fprintf(stderr, "input sizes mismatch (%s)\n", filename_m);
               for (int cx = 0; cx < nvertices; cx++)
               {
                       int i = mesh.v[cx].ij[0];
                       int j = mesh.v[cx].ij[1];
                       mesh.v[cx].im[ni].i = (m[2*(i+j*w)]  /wimax + ni)/nimages;
                       // abscisse dans l'atlas : située pour l'image ni dans [ni/nimages; (ni+1)/nimages]. On ajoute donc à ni/nimages, la coordonnée dans l'image (m[2*(i+j*w)]) multiplié par la part attribué à chaque pixel : (1/(wimax*nimages)).
                       mesh.v[cx].im[ni].j = m[2*(i+j*w)+1]/himax;
               }
               free(m);
       }

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
                                .v1 = q[1], .v2 = q[0]};
                        mesh.f[cx+1] = (struct face) {.v0 = q[3], 
                                .v1 = q[2], .v2 = q[1]};
			cx += 2;
		}
	}
	assert(cx == nfaces);

        // pour chaque face, calculer la normale et choisir la meilleure image
        for (int i = 0; i < nfaces; i++)
        {
                mesh.f[i].isvisible = malloc(mesh.nimages * sizeof(bool));
                struct face mf = mesh.f[i];
                triangle_normal(mesh.f[i].n, mesh.v[mf.v0].xyz, mesh.v[mf.v1].xyz, mesh.v[mf.v2].xyz);
                double sp = -1;
                double c_n[3];
                for (int ni = 0; ni < nimages; ni++)
                {
                        mesh.f[i].isvisible[ni] = ith_face_is_visible_in_image(mesh, i, ni);
                        if (ith_face_is_visible_in_image(mesh, i, ni))
                        {
                                for (int l = 0; l < 3; l++)
                                        c_n[l] = cam_n[3*ni+l];
                                if (fabs(scalar_product(c_n, mesh.f[i].n, 3)) > 1)
                                        printf("wARNING: scalar product error\n");
                                if (fabs(scalar_product(c_n, mesh.f[i].n, 3)) > sp)
                                {
                                        sp = fabs(scalar_product(c_n, mesh.f[i].n, 3));
                                        mesh.f[i].im = ni;
                                }
                        }
                }
                printf("face %d image %d\n", i, mesh.f[i].im);
        }
        // problème, on choisit la vue qui est la plus en face de la face. Cependant cette face peut être cachée sur cette vue. Il faut donc rafiner le choix et commencer par déterminer sur quelles vues la face est visible.

        write_ply_t(filename_ply, filename_a, mesh);
        write_ply_map_t(filename_ply, filename_a, mesh);
	return 0;
}

int main(int c, char *v[]) { return main_colormultiple(c,v); }

