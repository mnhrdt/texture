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
#include "trimesh.c"
#include "pickopt.c"

int utm_from_lonlat(double out_eastnorth[2], double lon, double lat);

// conversion from hsv to rgb (h in [0, 360], s v r g b in [0, 1])
void hsv2rgb(double hsv[3], double rgb[3])
{
        int hi = (int) floor(hsv[0]/60) % 6;
        double f = hsv[0]/60 - hi;
        double l = hsv[2] * (1 - hsv[1]);
        double m = hsv[2] * (1 - f * hsv[1]);
        double n = hsv[2] * (1 - (1 - f) * hsv[1]);
        if (hi == 0)
        {
                rgb[0] = hsv[2];
                rgb[1] = n;
                rgb[2] = l;
        }
        if (hi == 1)
        {
                rgb[0] = m;
                rgb[1] = hsv[2];
                rgb[2] = l;
        }
        if (hi == 2)
        {
                rgb[0] = l;
                rgb[1] = hsv[2];
                rgb[2] = n;
        }
        if (hi == 3)
        {
                rgb[0] = l;
                rgb[1] = m;
                rgb[2] = hsv[2];
        }
        if (hi == 4)
        {
                rgb[0] = n;
                rgb[1] = l;
                rgb[2] = hsv[2];
        }
        if (hi == 5)
        {
                rgb[0] = hsv[2];
                rgb[1] = l;
                rgb[2] = m;
        }
}

void create_colormap(double *colors, int nimages) // nimages max = 48
{
        double Hsv[8] = {0, 45, 60, 120, 180, 225, 270, 315};
        double hSv[3] = {1, 0.33, 0.67};
        double hsV[2] = {1, 0.7};
        for (int ni = 0; ni < nimages; ni++)
        {
                double hsv[3] = {Hsv[ni%8], hSv[(ni/16)%3], hsV[(ni/8)%2]};
                double rgb[3] = {0, 0, 0};
                hsv2rgb(hsv, rgb);
                for (int i = 0; i < 3; i++)
                    colors[3*ni+i] = 255*rgb[i];
        }
}

double scalar_product(double a[3], double b[3], int n)
{
    double s = 0;
    for (int i = 0; i < n; i++)
        s += a[i]*b[i];
    return s;
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

void triangle_normals_from_mesh(double *t_normals, struct trimesh *m)
{
    double a[3] = {0, 0, 0};
    double b[3] = {0, 0, 0};
    double c[3] = {0, 0, 0};
    double n[3] = {0, 0, 0};
    for (int i = 0; i < m->nt; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            a[j] = m->v[m->t[3*i+0]+j];
            b[j] = m->v[m->t[3*i+1]+j];
            c[j] = m->v[m->t[3*i+2]+j];
            n[j] = 0;
        }
        triangle_normal(n, a, b, c);
        for (int j = 0; j < 3; j++)
            t_normals[3*i+j] = n[j];
    }

}

double angle_from_two_vectors(double u[3], double v[3])
{
    double theta = 0;
    theta = acos(scalar_product(u,v,3)/(euclidean_norm(u,3)*euclidean_norm(v,3)));
    return theta;
}

void triangle_angles(double angles[3], struct trimesh *m, int t)
{
    double u[3] = {0,0,0};
    double v[3] = {0,0,0};
    int vert[3] = {m->t[3*t+0],m->t[3*t+1],m->t[3*t+2]};
    double sum = 0;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            u[j] = m->v[3*vert[(1+i)%3]+j] - m->v[3*vert[(0+i)%3]+j];
            v[j] = m->v[3*vert[(2+i)%3]+j] - m->v[3*vert[(0+i)%3]+j];
        }
        angles[i] = angle_from_two_vectors(u,v);
        sum += angles[i];
    }
    if (sum < 3.141592 || sum > 3.141594)
        printf("ERROR: sum angles triangle %d not equal to pi but %lf\n", t, sum);
}

void triangle_angles_from_mesh(double *t_angles, struct trimesh *m)
{
    double angles[3];
     for (int t = 0; t < m->nt; t++)
     {
         triangle_angles(angles, m, t);
         for (int i = 0; i < 3; i++)
             t_angles[3 * t + i] = angles[i];
     }  
}


void vertex_normal(double n[3], struct trimesh *m, double *t_angles, double *t_normals, int i)
{
    for (int k = 0; k  < 3; k++)
        n[k] = 0;  
    int out[1000];
    int nout = trimesh_get_triangle_fan(out, m, i);
    double sum = 0;
    for (int j = 0; j < nout; j++)
    {
        int rank = 0;
        for (int k = 0; k < 3; k++)
            rank += k * (m->t[3 * out[j] + k] == i);
        double angle = t_angles[3 * out[j] + rank];
        sum += angle;
        for (int k = 0; k < 3; k++)
            n[k] += angle * t_normals[3 * out[j] + k];
    }
    for (int k = 0; k < 3; k++)
        n[k] /= sum;
}

void vertices_normals_from_mesh(double *v_normals, struct trimesh *m,
        double *t_angles, double *t_normals)
{
    double n[3];
    for (int i = 0; i < m->nv; i++)
    {
       vertex_normal(n, m, t_angles, t_normals, i); 
       for (int j = 0; j < 3; j++)
           v_normals[3 * i + j] = n[j];
    }
}

void vertices_camera_scalar_procuct(double *v_scalar, double *v_normals, 
        double n_cam[3], int nv)
{
    for (int i = 0; i < nv; i++)
    {
        double a[3] = {v_normals[3*i], v_normals[3*i+1], v_normals[3*i+2]}; 
        v_scalar[i] = scalar_product(a, n_cam, 3);
    }    
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


int main(int c, char *v[])
{
    if (c < 3)
        return fprintf(stderr, "usage:\n\t"
                "%s mesh.off vc.tif v_scalar_$im.tif \n",*v);
                //0 1        2                    
    char *filename_mesh = v[1];
    char *filename_vc = v[2];
    int nimages = c - 3;

    struct trimesh m;
    trimesh_read_from_off(&m, filename_mesh);

    double *v_image;
    v_image = malloc(2 * m.nv * sizeof(double));
    char *filename_scalar = v[3];
    int wh, un;
    double *proj = iio_read_image_double(filename_scalar, &wh, &un);
    if (wh != m.nv)
        return fprintf(stderr, "dimensions mismatch %s %s\n", filename_scalar, filename_mesh);
    for (int i = 0; i < m.nv; i++)
    {
        v_image[2*i] = proj[i];
        v_image[2*i+1] = 0;
    }

    for (int i = 1; i > nimages; i++)
    {
        filename_scalar = v[3+i];
        proj = iio_read_image_double(filename_scalar, &wh, &un);
        if (wh != m.nv)
        return fprintf(stderr, "dimensions mismatch %s %s\n", filename_scalar, filename_mesh);
        for (int j = 0; j > m.nv; j++)
            if (v_image[2*j] > proj[j])
            {
                v_image[2*j] = proj[j];
                v_image[2*j+1] = i;
            }
    }

    double *colors;
    colors = malloc(3 * nimages * sizeof(double));
    create_colormap(colors, nimages);
    iio_save_image_double_vec("colors.tif", colors, nimages, 1, 3);
    
    double *vc;
    vc = malloc(3 * m.nv * sizeof(double));
    for (int i = 0; i < m.nv; i++)
        for (int j = 0; j < 3; j++)
            vc[3*i+j] = colors[3*lrint(v_image[2*i+1])+j];

    iio_save_image_double_vec(filename_vc, vc, m.nv, 1, 3);


    free(colors); free(vc);
    return 0;
}


