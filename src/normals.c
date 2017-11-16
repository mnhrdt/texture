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

int main(int argc, char *v[])
{
    if (argc < 3)
        return fprintf(stderr, "usage:\n\t"
                "%s mesh.off triangles_normal.tif\n",*v);
                //0 1        2                    
    char *filename_mesh = v[1];
    char *filename_n = v[2];
    char *filename_a = v[3];


    struct trimesh m;
    trimesh_read_from_off(&m, filename_mesh);
    trimesh_fill_triangle_fans(&m);

    double *t_normals;
    t_normals = malloc(3 * m.nt * sizeof(double));
    triangle_normals_from_mesh(t_normals, &m);
    printf("normales calculées\n");
    
    double *t_angles;
    t_angles = malloc(3 * m.nt * sizeof(double));
    triangle_angles_from_mesh(t_angles, &m);
    printf("angles calculées\n");

    double *v_normals;
    v_normals = malloc(3 * m.nv * sizeof(double));
    vertices_normals_from_mesh(v_normals, &m, t_angles, t_normals);
    printf("assignation des angles aux sommets\n");


    free(t_normals); free(t_angles); free(v_normals);
    return 0;
}


