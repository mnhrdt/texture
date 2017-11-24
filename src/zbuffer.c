#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

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

void fill_vector_from_file(double *xywh, char *filename, int l)
{
    // read the informations about the crop
    FILE *file;
    file = fopen(filename,"r");
    if (!file)
        fprintf(stderr, "fopen(%s) failed\n", filename);
        

    for (int i = 0; i < l; i++)
        if ((fscanf(file, "%lf", &xywh[i])) != 1)
             fprintf(stderr, "could not read element %d of %s\n", i, filename);
    fclose(file);
}


bool fill_three_point_of_img_copy(struct trimesh *m,
        double origin[3],
        int signed_zone,
        struct rpc *huge_rpc,
        double xywh[4],
        double *img_copy,
        bool *v_visibility,
        int v_coord_im[3][2],
        int nt)
{
    int wi = lrint(xywh[2]);
    int hi = lrint(xywh[3]);
    
    int vertices[3] = {m->t[3*nt+0], m->t[3*nt+1], m->t[3*nt+2]};
    bool in_image = true;
    for (int l = 0; l < 3; l++)
    {
        // get vertices coordinates in pixel and utm
        double e = m->v[3 * vertices[l] + 0] + origin[0];
        double n = m->v[3 * vertices[l] + 1] + origin[1];
        double z = m->v[3 * vertices[l] + 2] + origin[2];
        if (z<-100)
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
            ij_int[k] = (int) floor(ij[k]) - floor(xywh[k]);
        if (!is_in_crop_int(ij_int, xywihi)) 
        {   
            //    printf("ij_int : %d %d\n", ij_int[0], ij_int[1]);
            in_image = false; 
            return in_image;
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
            img_copy[3*(ij_int[1]*wi+ij_int[0])+2] = z-origin[2];
        }
    }
    return in_image;

}

void fill_img_copy(struct trimesh *m,
        double origin[3],
        int signed_zone, 
        struct rpc *huge_rpc,
        double xywh[4],
        double *img_copy,
        bool *v_visibility,
        double n_cam[3])
{
    int wi = lrint(xywh[2]);
    int hi = lrint(xywh[3]);
    for (int nt = 0; nt < m->nt; nt++)
    {
        int vertices[3] = {m->t[3*nt+0], m->t[3*nt+1], m->t[3*nt+2]};
        double v_coord_scaled[3][3];
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
                v_coord_scaled[i][j] = (m->v[3 * vertices[i] + j]);
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

        int v_coord_im[3][2];
        // loop over the vertices of the well-oriented triangle.
        bool in_image = fill_three_point_of_img_copy(m, origin, signed_zone, huge_rpc, xywh, img_copy, v_visibility, v_coord_im, nt);

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
}

void check_vertex_visibility_and_fill_output(struct trimesh *m, 
        double origin[3], 
        int signed_zone, 
        struct rpc *huge_rpc, 
        double xywh[4], 
        double *img_copy, 
        double *out, 
        bool *v_visibility, 
        double n_cam[3], 
        double resolution)
{
    int wi = lrint(xywh[2]);
    // loop over triangle vertices and give them coordinate of projection on 
    // image only if there are no occlusions (they are part of a well oriented
    // triangle and are close enough to the highest point projected on 
    // the image.
    for (int nv = 0; nv < m->nv; nv++)
        // check if the vertex can be projected on the image
        // and is part of a well oriented triangle.
        // Otherwise, go to the next vertex.
    {

        if (!v_visibility[nv])
        {
            for (int i = 0; i < 3; i++)
                out[3*nv+i] = NAN;
            continue;
        }  

        double e = m->v[3*nv+0] + origin[0];
        double n = m->v[3*nv+1] + origin[1];
        double z = m->v[3*nv + 2]+origin[2];
        if (z<-100)
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
        {
            out[3*nv+0] = lonlat[0]; 
            out[3*nv+1] = lonlat[1]; 
            out[3*nv+2] = z; 
        }
        else
        {
            //printf("img_copy %lf z %lf\n", img_copy[3*(ij_int[1]*wi + ij_int[0]) + 2], z);
            double enz[3] = {e-origin[0], n-origin[1], z};
            double enz_ref[3] = {0, 0, 0};
            for (int k = 0; k < 3; k++)
                enz_ref[k] = img_copy[3*(ij_int[1]*wi+ij_int[0])+k];
            // Compute squared distance to saved point in img_copy.
            // Substract squared projection on camera direction.
            double diff = 0;
            for (int k = 0; k < 3; k++)
                //diff += pow(enz[k]-enz_ref[k],2)
                 //   * (1-pow(n_cam[k],2));
                diff += pow(enz[k]-enz_ref[k],2) - pow((enz[k]-enz_ref[k])*n_cam[k],2);
            // By Pythagore the result should be less than 0.3^2
            // (resolution is 30cm per pixel
 //           printf("diff %lf enz %lf %lf %lf t %lf\n", diff, enz[0], enz[1], enz[2], 2 * pow(resolution,2));
//            printf("img_copy %lf %lf %lf\n", enz_ref[0], enz_ref[1], enz_ref[2]); 


            if (diff < 2 * pow(resolution,2))
            {
                out[3*nv+0] = lonlat[0]; 
                out[3*nv+1] = lonlat[1]; 
                out[3*nv+2] = z; 
            }

        }
    }
}

void scale_and_origin_from_file(double scale[3], double origin[3], char *filename_scale)
{
    double xywh[4];
    fill_vector_from_file(xywh, filename_scale, 4);
    for (int i = 0; i < 2; i++)
    {
        scale[i] = xywh[i];
        origin[i] = xywh[2+i];
    }
}

int main_zbuffer_meter(int c, char *v[])
{
    double resolution = atof(pick_option(&c, &v, "-res", "0.3"));
if (c < 7)
        return fprintf(stderr, "usage:\n\t"
                "%s scale.txt zone rpc_i xywh_i.txt mesh.off out.tif\n", *v);
                //0 1         2    3     4          5        6       
    char *filename_scale = v[1];
    int signed_zone = atoi(v[2]);
    char *filename_rpc = v[3];
    char *filename_corners = v[4];
    char *filename_mesh = v[5];
    char *filename_out = v[6];



    // read georeferencing transform data
    double origin[3] = {0, 0, 0};
    double scale[3] = {1, 1, 1};
    scale_and_origin_from_file(scale, origin, filename_scale);

    struct trimesh m[1];
    trimesh_read_from_off(m, filename_mesh);

    // read the input image (small jpg, png or tif crop corresponding to the DSM)
//    int wi, hi;
//    float *img = iio_read_image_float(filename_img, &wi, &hi);
//    if (!img)
//        return fprintf(stderr, "iio_read(%s) failed\n", filename_img);

    // read the informations about the crop
    double xywh[4];
    fill_vector_from_file(xywh, filename_corners, 4);

    int wi = lrint(xywh[2]);
    int hi = lrint(xywh[3]);
    // allocate space for same size multispectral image for (e, n, z)
    // initialise height to -100.
    double *img_copy = malloc(3 * wi * hi * sizeof(double));
    for (int i=0; i < wi*hi; i++)
    {
        for (int j = 0; j < 2; j++)
            img_copy[3*i+j] = 0;
        img_copy[3*i+2] = -1000;
    }
    
    // allocate space for output image and initialise to NAN.
    double *out;
    out = malloc(3 * m->nv * sizeof(double));
    for (int i = 0; i < 3 * m->nv; i++)
        out[i] = NAN;

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
    fill_img_copy(m, origin, signed_zone, huge_rpc, xywh, img_copy, 
            v_visibility, n_cam);
    check_vertex_visibility_and_fill_output(m, origin, signed_zone, huge_rpc, 
            xywh, img_copy, out, v_visibility, n_cam, resolution);

    iio_save_image_double_vec("exp/soutput/image_copy.tif", img_copy, wi, hi, 3);
    iio_save_image_double_vec(filename_out, out, m->nv, 1, 3);
    free(img_copy); free(out); free(v_visibility);
    return 0;
}
int main(int c, char *v[]) { return main_zbuffer_meter(c,v); }
