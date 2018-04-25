#ifndef NORMALS_H
#define NORMALS_H

#include "rpc.h"
#include "trimesh.h"

int utm_from_lonlat(double out_eastnorth[2], double lon, double lat);

double scalar_product(double a[3], double b[3], int n);

void cross_product(double axb[3], double a[3], double b[3]);

double euclidean_norm(double *x, int n);

// coordinates of the normal to a triangle in a 3D space
void triangle_normal(double n[3], double a[3], double b[3], double c[3]);
    // les sommets sont donnés dans le sens direct

void triangle_normals_from_mesh(double *t_normals, struct trimesh *m);

double angle_from_two_vectors(double u[3], double v[3]);

void triangle_angles(double angles[3], struct trimesh *m, int t);

void triangle_angles_from_mesh(double *t_angles, struct trimesh *m);


void vertex_normal(double n[3], struct trimesh *m, double *t_angles, double *t_normals, int i);

void vertices_normals_from_mesh(double *v_normals, struct trimesh *m,
        double *t_angles, double *t_normals);

void vertices_camera_scalar_product(double *v_scalar, double *v_normals, 
        double n_cam[3], int nv);


#endif

