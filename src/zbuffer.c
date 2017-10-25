#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gdal.h>

#include "iio.h"
#include "drawtriangle.c"
#include "trimesh.c"
#include "rpc.c"
#include "pickopt.c"

int utm_from_lonlat(double out_eastnorth[2], double lon, double lat);
void lonlat_from_eastnorthzone(double out_lonlat[2], double e, double n, int z);

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

// coordinates of the normal to a triangle in a 3D space
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

// get the satellite direction using only the rpc data
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
    double norm = euclidean_norm(n, 3);
    for (int i = 0; i < 3; i++) n[i] /= norm;

    norm = euclidean_norm(n, 3);
    if (norm < 0.9999 || norm > 1.0000001)
        printf("WARNING: normalisation error in camera_direction, norme = %.16lf\n", norm);
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
    double *img_copy;
    double c[3][3]; // vertices utm coordinates (without origin) !!!! PROBLEM ?
};

static bool is_large_triangle(float abc[3][2])
{
    bool a = false;
    for (int i = 0; i < 3; i++)
        a |= ((abc[i][0] - abc[(i+1)%3][0])>1 || (abc[i][1] - abc[(i+1)%3][1])>1);

    a &= (abc[0][1] != abc[1][1] || abc[0][1] != abc[2][1] || abc[1][1] != abc[2][1]);
    return a;
}

static bool is_in_crop_int(int *ij, int *xywh)
{
    return ((ij[0]>= xywh[0]) && (ij[0] < xywh[0] + xywh[2]) 
            && (ij[1]>= xywh[1]) && (ij[1] < xywh[1] + xywh[3])); 
}

static bool is_in_crop(double *ij, double *xywh)
{
    return ((ij[0]>= xywh[0]) && (ij[0] <= xywh[0] + xywh[2]) 
            && (ij[1]>= xywh[1]) && (ij[1] <= xywh[1] + xywh[3])); 
}

static void interpolate_vertices_values(int i, int j, void *ee)
{
    struct im_ver *e = ee;
    int w = e->w;
    int h = e->h;
    double c[3][3];
    for (int k = 0; k < 3; k++)
        for (int l = 0; l < 3; l++)
            c[k][l] = e->c[k][l];
    int u[2] = {e->v[2][0] - e->v[0][0], e->v[2][1] - e->v[0][1]};
    int v[2] = {e->v[1][0] - e->v[0][0], e->v[1][1] - e->v[0][1]};
    // vérifie que l'on a un vrai triangle abc
    // en cas de triangle dégénéré, une droite,
    if ((u[1]*v[0] - u[0]*v[1])==0)
    {
        // si b=a ou b=c 
        if ((e->v[0][0] == e->v[1][0] && e->v[0][1] == e->v[1][1]) ||
                (e->v[2][0] == e->v[1][0] && e->v[2][1] == e->v[1][1]    )) 
        {
            // on vérifie que a!=c
            assert(e->v[0][0] != e->v[2][0] || e->v[0][1] != e->v[2][1]);
            // puis on remplit la droite
            double ax[2] = {i-e->v[0][0], j - e->v[0][1]};
            double ac[2] = {e->v[2][0] - e->v[0][0], e->v[2][1] - e->v[0]    [1]};
            // en faisant à chaque point le test de hauteur.
            // On ne remplit l'image que si le nouveau point est plus 
            // haut.
            if (c[0][2] + (c[2][2] - c[0][2])*euclidean_norm(ax, 2) 
                    / euclidean_norm(ac, 2) 
                    > e->img_copy[3 * (i+j*w) + 2])
                // pour chaque coordonnée
                for (int k = 0; k < 3; k++)
                    e->img_copy[3*(i+j*w) + k] = c[0][k] 
                        + (c[2][k] - c[0][k])*euclidean_norm(ax, 2) 
                        / euclidean_norm(ac, 2);
        }
        // si a=c
        if (e->v[0][0] == e->v[2][0] && e->v[0][1] == e->v[2][1])
        {
            // on vérifie que b!=a
            assert(e->v[0][0] != e->v[1][0] || e->v[0][1] != e->v[1][1]);
            double ax[2] = {i-e->v[0][0], j - e->v[0][1]};
            double ac[2] = {e->v[1][0] - e->v[0][0], e->v[1][1] - e->v[0]    [1]};
            // et on remplit la droite avec test
            if (c[0][2] + (c[1][2] - c[0][2])*euclidean_norm(ax, 2)  
                    / euclidean_norm(ac, 2) 
                    > e->img_copy[3 * (i+j*w) + 2])
                // pour chaque coordonnée
                for (int k = 0; k < 3; k++) 
                    e->img_copy[3*(i+j*w)+k] = c[0][k] 
                        + (c[1][k] - c[0][k])*euclidean_norm(ax, 2) 
                        / euclidean_norm(ac, 2);

        }
    }
    // si triangle normal, on interpole de façon classique
    else
    {
        if (i>=w)
            printf("w : %d, i : %d\n", w, i);
        assert(i<w);
        assert(j<h);
        assert(3*(i+j*w)+2 < 3*w*h);
        double alpha = (double) ((j-e->v[0][1])*v[0] - (i-e->v[0][0])*v[1])
            / (u[1]*v[0] - u[0]*v[1]);
        double beta = (double) (-(j-e->v[0][1])*u[0] + (i-e->v[0][0])*u[1])
            / (u[1]*v[0] - u[0]*v[1]);
        // test pour savoir si nouveau point plus haut
        if (c[0][2] + alpha * (c[2][2] - c[0][2]) 
                + beta * (c[1][2] - c[0][2]) 
                > e->img_copy[3*(i+j*w) + 2])
            // remplissage de chaque coordonnée
            for (int k = 0; k < 3; k++)
                e->img_copy[3*(i+j*w) + k] = c[0][k] 
                    + alpha * (c[2][k] - c[0][k]) 
                    + beta * (c[1][k] - c[0][k]);
    }

}


static void get_scale_and_origin_from_gdal(double scale[2], double origin[2],
        char *filename_dsm)
{
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
}


int main_zbuffer(int c, char *v[])
{
    if (c < 8)
        return fprintf(stderr, "usage:\n\t"
                "%s dsm.tif zone img_i.tif rpc_i xywh_i.txt mesh.off out.tif\n", *v);
                //0 1       2    3         4     5          6       7
    double resolution = atof(pick_option(&c, &v, "-res", "0.3"));
    char *filename_dsm = v[1];
    int signed_zone    = atoi(v[2]);
    char *filename_corners = v[5];
    char *filename_img = v[3];
    char *filename_rpc= v[4];
    char *filename_out = v[7];
    char *filename_mesh = v[6];

    // read the whole input DSM (typically, rather small)
    int w, h;
    float *x = iio_read_image_float(filename_dsm, &w, &h);
    if (!x)
        return fprintf(stderr, "iio_read(%s) failed\n", filename_dsm);

    // read georeferencing transform using GDAL
    double origin[2] = {0, 0};
    double scale[2] = {1, 1};
    get_scale_and_origin_from_gdal(scale, origin, filename_dsm);

    struct trimesh m[1];
    trimesh_read_from_off(m, filename_mesh);

    // read the input image (small jpg, png or tif crop corresponding to the DSM)
    int wi, hi;
    float *img = iio_read_image_float(filename_img, &wi, &hi);
    if (!img)
        return fprintf(stderr, "iio_read(%s) failed\n", filename_img);

    // read the informations about the crop
    FILE *corners;
    corners = fopen(filename_corners,"r");
    if (!corners)
        return fprintf(stderr, "fopen(%s) failed\n", filename_corners);

    double xywh[4];
    for (int i = 0; i < 4; i++)
        if ((fscanf(corners, "%lf", &xywh[i])) != 1)
            return fprintf(stderr, "could not read element %d of %s\n", i, filename_corners);
    fclose(corners);


    // allocate space for same size multispectral image for (e, n, z)
    // initialise height to -100.
    double *img_copy = malloc(3 * wi * hi * sizeof(double));
    for (int i=0; i < wi*hi; i++)
    {
        for (int j = 0; j < 2; j++)
            img_copy[3*i+j] = 0;
        img_copy[3*i+2] = -100;
    }

    // allocate space for output image and initialise to -1.
    //        double *out = malloc(2 * w * h * sizeof(double));
    //        for (int i = 0; i < 2 * w * h; i++)
    //                out[i] = -1;
    double *out = malloc(2 * m->nv * sizeof(double));
    for (int i = 0; i < 2 * m->nv; i++)
        out[i] = -1;

    // allocate space for vertices visibility. Initialise to false.
    bool *v_visibility = malloc(m->nv * sizeof(bool));
    for (int i = 0; i < m->nv; i++)
        v_visibility[i] = false;

    struct rpc huge_rpc[1];
    read_rpc_file_xml(huge_rpc, filename_rpc);

    // get camera direction
    double n_cam[3];
    camera_direction(n_cam, huge_rpc);


    // loop over mesh faces and fill img_copy with height of lidar point
    // keep only the highest heigtht for each point to deal with occlusions
    for (int nt = 0; nt < m->nt; nt++)
    {
        int vertices[3] = {m->t[3*nt+0], m->t[3*nt+1], m->t[3*nt+2]};
        double v_coord_scaled[3][3];
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 2; j++)
                v_coord_scaled[i][j] = m->v[3 * vertices[i] + j] * scale[j];
            v_coord_scaled[i][2] = m->v[3 * vertices[i] + 2];
        }

        // get triangle normal
        double n_triangle[3];
        triangle_normal(n_triangle, v_coord_scaled[0], 
                v_coord_scaled[1], v_coord_scaled[2]);

        // test if the triangle normal is in the opposite direction of the 
        // camera. If same direction, go to the next triangle.
        double tmp_d = 0;
        for (int i = 0; i < 3; i++)
            tmp_d += n_triangle[i] * n_cam[i];
        if (tmp_d < 0)
            continue;

        bool in_image = true;
        int v_coord_im[3][2];
        // loop over the vertices of the well-oriented triangle.
        for (int l = 0; l < 3; l++)
        {
            // get vertices coordinates in pixel and utm
            double i = m->v[3 * vertices[l] + 0];
            double j = m->v[3 * vertices[l] + 1];
            double e = i * scale[0] + origin[0];
            double n = j * scale[1] + origin[1];
            double z = m->v[3 * vertices[l] + 2];
            if (z<0)
                z = 20.0; // TO DO: put mean value

            // get projection coordinates in huge image
            double lonlat[2] = {0, 0};
            lonlat_from_eastnorthzone(lonlat, e, n, signed_zone);
            double lonlatheight[3] = {lonlat[0], lonlat[1], z};
            double ij[2] = {0, 0};
            rpc_projection(ij, huge_rpc, lonlatheight);

            // check if the vertex is in the cropped image
            // else put NaN in output and go to the next vertex.
            int ij_int[2];
            int xywihi[4] = {0, 0, wi, hi};
            for (int k = 0; k < 2; k++)
                // projection sur grande image,
                // il faut donc soustraire les coord du coin 
                // de l'image découpée.
                ij_int[k] = (int) round(ij[k]) - round(xywh[k]);
            if (!is_in_crop_int(ij_int, xywihi)) 
            {
                //    printf("ij_int : %d %d\n", ij_int[0], ij_int[1]);
                in_image = false;
                for (int k = 0; k < 2; k++)
                    out[2*vertices[l] + k] = NAN;
                //out[2*((int) round(i)+(int) round(j)*w) + k] = NAN;
                continue;
            }
            // indicate that the vertex is part of a well-oriented face
            assert(vertices[l] < m->nv);
            assert(vertices[l] >= 0);
            v_visibility[vertices[l]] = true;

            for (int k = 0; k < 2; k++)
                v_coord_im[l][k] = ij_int[k];

            // fill img_copy only if it is the highest point
            if (img_copy[3*(ij_int[1]*wi+ij_int[0])+2] <= z)
            {
                img_copy[3*(ij_int[1]*wi+ij_int[0])+0] = e-origin[0];
                img_copy[3*(ij_int[1]*wi+ij_int[0])+1] = n-origin[1];
                img_copy[3*(ij_int[1]*wi+ij_int[0])+2] = z;
            }
        }

        if (!in_image)
            continue;
        // fill in the triangle formed in img_copy.
        struct im_ver e = {.w = wi, .h = hi, .img_copy = img_copy};
        float abc[3][2];
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                abc[i][j] = (float) v_coord_im[i][j];
                e.v[i][j] = v_coord_im[i][j];
                e.c[i][j] = v_coord_scaled[i][j];
            }
            e.c[i][2] = v_coord_scaled[i][2];
        }
        //if (visible && is_large_triangle(abc))
        if (is_large_triangle(abc))
        {
            traverse_triangle(abc, interpolate_vertices_values, &e);
        }
    }
//    iio_save_image_double_vec("essai_curve/data/img_copy.tif", img_copy, wi, hi,3);


    // loop over triangle vertices and give them coordinate of projection on 
    // image only if there are no occlusions (they are part of a well oriented
    // triangle and are close enough to the highest point projected on 
    // the image.
    double enz_ref[3];
    double diff;
    for (int nv = 0; nv < m->nv; nv++)
        // check if the vertex can be projected on the image
        // and is part of a well oriented triangle.
        // Otherwise, go to the next vertex.
    {

        int i = (int) round(m->v[3*nv + 0]);
        int j = (int) round(m->v[3*nv + 1]);

        if (!v_visibility[nv])
            continue;

        double e = i*scale[0] + origin[0];
        double n = j*scale[1] + origin[1];
        double z = m->v[3*nv + 2];
        if (z<0)
            z = 20.0; // TO DO: put mean value

        // get projection coordinates in huge image
        double lonlat[2] = {0, 0};
        lonlat_from_eastnorthzone(lonlat, e, n, signed_zone);
        double lonlatheight[3] = {lonlat[0], lonlat[1], z};
        double ij[2] = {0, 0};
        rpc_projection(ij, huge_rpc, lonlatheight);

        // round the indices
        int ij_int[2];
        for (int k = 0; k < 2; k++)
            ij_int[k] = (int) round(ij[k]) - round(xywh[k]);
        if (img_copy[3*(ij_int[1]*wi + ij_int[0]) + 2] == z)
            for (int k = 0; k < 2; k++)
                out[2*nv+k] = ij[k]-xywh[k];
        //out[2*(i+j*w)+k] = ij[k]-xywh[k];
        else
        {
            //printf("img_copy %lf z %lf\n", img_copy[3*(ij_int[1]*wi + ij_int[0]) + 2], z);
            double enz[3] = {e-origin[0], n-origin[1], z};
            for (int k = 0; k < 3; k++)
                enz_ref[k] = img_copy[3*(ij_int[1]*wi+ij_int[0])+k];
            // Compute squared distance to saved point in img_copy.
            // Substract squared projection on camera direction.
            diff = 0;
            for (int k = 0; k < 3; k++)
                diff += pow(enz[k]-enz_ref[k],2)
                    * (1-pow(n_cam[k],2));
            // By Pythagore the result should be less than 0.3^2
            // (resolution is 30cm per pixel
            if (diff < 2 * pow(resolution,2))
                for (int k = 0; k < 2; k++)
                    out[2*nv+k] = ij[k] - xywh[k];
            //out[2*(i+j*w)+k] = ij[k] - xywh[k];

        }
    }
    iio_save_image_double(filename_out, out, m->nv, 2);
    //iio_save_image_double_vec(filename_out, out, w, h, 2);

    free(img_copy); free(out); free(v_visibility);
    return 0;
}

int main(int c, char *v[]) { return main_zbuffer(c,v); }
