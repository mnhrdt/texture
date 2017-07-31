#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gdal.h>

#include "iio.h"
#include "drawtriangle.c"

double euclidean_norm(double *x, int n)
{
        return n > 0 ? hypot(*x, euclidean_norm(x+1, n-1)) : 0;
}

// extrapolate by nearest neighbour
static void keep_in_bounds(int w, int h, int *ij)
{
                if (ij[0] < 0) ij[0] = 0;
                if (ij[1] < 0) ij[1] = 0;
                if (ij[0] >= w) ij[0] = w-1;
                if (ij[1] >= h) ij[1] = h-1;
}

struct im_ver{
        int w, h;    // lidar width, height
        int v[3][2]; // coordonnées i j des sommets dans image
        float *img_copy;
        float z[3]; // valeur des sommets
};

static bool is_large_triangle(float abc[3][2])
{
        bool a = false;
        for (int i = 0; i < 3; i++)
                a |= ((abc[i][0] - abc[(i+1)%3][0])>1 || (abc[i][1] - abc[(i+1)%3][1])>1);

        a &= (abc[0][1] != abc[1][1] || abc[0][1] != abc[2][1] || abc[1][1] != abc[2][1]);
        return a;
}

static void interpolate_vertices_values(int i, int j, void *ee)
{
        struct im_ver *e = ee;
        int w = e->w;
        int h = e->h;
        float z[3];
        for (int l = 0; l < 3; l++)
                z[l] = e->z[l];
        int u[2] = {e->v[2][0] - e->v[0][0], e->v[2][1] - e->v[0][1]};
        int v[2] = {e->v[1][0] - e->v[0][0], e->v[1][1] - e->v[0][1]};
        // vérifie que l'on a un vrai triangle
        // en cas de triangle dégénéré
        if ((u[1]*v[0] - u[0]*v[1])==0)
        {
                if ((e->v[0][0] == e->v[1][0] && e->v[0][1] == e->v[1][1]) ||
                                (e->v[2][0] == e->v[1][0] && e->v[2][1] == e->v[1][1]    )) 
                {
                        assert(e->v[0][0] != e->v[2][0] || e->v[0][1] != e->v[2][1]);
                        double ax[2] = {i-e->v[0][0], j - e->v[0][1]};
                        double ac[2] = {e->v[2][0] - e->v[0][0], e->v[2][1] - e->v[0]    [1]};
                        if (z[0] + (z[2] - z[0])*euclidean_norm(ax,     2) / euclidean_norm(ac, 2) > e->img_copy[i+j*w])
                                e->img_copy[i+j*w] = z[0] + (z[2] - z[0])*euclidean_norm(ax, 2) / euclidean_norm(ac, 2);
                }
                if (e->v[0][0] == e->v[2][0] && e->v[0][1] == e->v[2][1])
                {
                        assert(e->v[0][0] != e->v[1][0] || e->v[0][1] != e->v[1][1]);
                        double ax[2] = {i-e->v[0][0], j - e->v[0][1]};
                        double ac[2] = {e->v[1][0] - e->v[0][0], e->v[1][1] - e->v[0]    [1]};
                        if (z[0] + (z[1] - z[0])*euclidean_norm(ax,     2)  / euclidean_norm(ac, 2) > e->img_copy[i+j*w])
                                e->img_copy[i+j*w] = z[0] + (z[1] - z[0])*euclidean_norm(ax, 2)  / euclidean_norm(ac, 2);
                }
        }
        else
        {
                double alpha = (double) ((j-e->v[0][1])*v[0] - (i-e->v[0][0])*v[1])
                        / (u[1]*v[0] - u[0]*v[1]);
                double beta = (double) (-(j-e->v[0][1])*u[0] + (i-e->v[0][0])*u[1])
                        / (u[1]*v[0] - u[0]*v[1]);
                if (z[0] + alpha * (z[2] - z[0]) + beta * (z[1] - z[    0]) > e->img_copy[i+j*w])
                        e->img_copy[i+j*w] = z[0] + alpha * (z[2] - z[0]) + beta * (z[1] - z[0]);
        }

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

void initialize_mesh_from_lidar(struct mesh_t *mesh, char *filename_dsm)
{
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


        // read the whole input DSM (typically, rather small)
        int w, h;
        float *x = iio_read_image_float(filename_dsm, &w, &h);
        if (!x)
                fprintf(stderr, "iio_read(%s) failed\n", filename_dsm);



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
        mesh->nv = nvertices;

        // count number of valid square faces
        int nfaces = 0;
        for (int j = 0; j < h-1; j++)
        for (int i = 0; i < w-1; i++)
        {
                int q[4] = {vid[j][i], vid[j+1][i], vid[j+1][i+1], vid[j][i+1]};
                if (q[0] >= 0 && q[1] >= 0 && q[2] >= 0 && q[3] >= 0)
                        nfaces += 2;
        }
        mesh->nf = nfaces;

        // initialize vertices with their coordinates in space and in the lidar
        mesh->v = malloc(mesh->nv*sizeof(struct vertex));
        int cx = 0;
        for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++)
        {
                if (!isfinite(height[j][i])) continue;
                // compute lonlat from eastnorth = {p[0], p[1]}
                double e = i* scale[0];// + origin[0]; // easting
                double n = j* scale[1];// + origin[1]; // northing
                double z = height[j][i];             // height
                mesh->v[cx].xyz[0] = e;
                mesh->v[cx].xyz[1] = n;
                mesh->v[cx].xyz[2] = z;

                mesh->v[cx].ij[0] = i;
                mesh->v[cx].ij[1] = j;

                cx += 1;
        }
        assert(cx == nvertices);
        // for each square in lidar, create two triangular faces.
        // Vertices order is such that the face is oriented towards the outside.
        mesh->f = malloc(mesh->nf*sizeof(struct face));
        cx = 0;
        for (int j = 0; j < h-1; j++)
                for (int i = 0; i < w-1; i++)
                {
                        int q[4] = {vid[j][i], vid[j+1][i], vid[j+1][i+1], vid[j][i+1]};
                        if (q[0] >= 0 && q[1] >= 0 && q[2] >= 0 && q[3] >= 0)
                        {
                        if (fabs(x[i+j*w]-x[i+1+(j+1)*w]) >= fabs(x[i+1+j*w]-x[i+(j+1)*w]))
                        {
                                mesh->f[cx] = (struct face) {.v0 = q[3],
                                        .v1 = q[0], .v2 = q[1]};
                                mesh->f[cx+1] = (struct face) {.v0 = q[3],
                                        .v1 = q[1], .v2 = q[2]};
                                cx += 2;
                        }
                        else
                        {
                                mesh->f[cx] = (struct face) {.v0 = q[0],
                                        .v1 = q[1], .v2 = q[2]};
                                mesh->f[cx+1] = (struct face) {.v0 = q[0],
                                        .v1 = q[2], .v2 = q[3]};
                                cx += 2;
                        }
                        }

                }
        assert(cx == nfaces);
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


int main_zbuffer(int c, char *v[])
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

        // create and initialize mesh with lidar
        struct mesh_t mesh;
        initialize_mesh_from_lidar(&mesh, filename_dsm);

        // loop over mesh faces and fill img_copy with height of lidar point
        // keep only the highest heigtht for each point to deal with occlusions
        for (int nf = 0; nf < mesh.nf; nf++)
        {
                bool visible = true;
                struct face mf = mesh.f[nf];
                int vertices[3] = {mf.v0, mf.v1, mf.v2};
                int v_coord[3][2];
                float v_z[3];
                for (int l = 0; l < 3; l++)
                {
                        int i = mesh.v[vertices[l]].ij[0];
                        int j = mesh.v[vertices[l]].ij[1];
                        double z = x[i+j*w];
                        v_z[l] = z;
                        if (z<0)
                                z = 20.0;
                        double xyz1[4] = {i, j, z, 1};
                        double ij_approx[3];
                        
                        matrix_product_mxn_nxp(ij_approx, P, xyz1, 3, 4, 1); 

                        int ij_int[2];
                        for (int k = 0; k < 2; k++)
                                ij_int[k] = (int) round(ij_approx[k]);

                        keep_in_bounds(wi, hi, ij_int);
                        v_coord[l][0] = ij_int[0];
                        v_coord[l][1] = ij_int[1];
                        if (img_copy[ij_int[1]*wi+ij_int[0]] <= z)
                                img_copy[ij_int[1]*wi+ij_int[0]] = z;
                        else
                                visible = false;
                }
                struct im_ver e = {.w = wi, .h = hi, .img_copy = img_copy};
                float abc[3][2];
                for (int i = 0; i < 3; i++)
                for (int j = 0; j < 2; j++)
                {
                        abc[i][j] = (float) v_coord[i][j];
                        e.v[i][j] = v_coord[i][j];
                        e.z[i] = v_z[i];
                }
                //if (visible && is_large_triangle(abc))
                if (is_large_triangle(abc))
                {
                        traverse_triangle(abc, interpolate_vertices_values, &e);
                }
        }
        iio_save_image_float("essai_curve/data/img_copy.tif", img_copy, wi, hi);


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

int main(int c, char *v[]) { return main_zbuffer(c,v); }
