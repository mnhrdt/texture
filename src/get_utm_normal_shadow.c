#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "iio.h"
#include "normals.h"
#include "drawtriangle.h"
#include "trimesh.h"
#include "rpc.h"
#include "pickopt.h"
#include "get_utm_normal_shadow.h"

int utm_from_lonlat(double out_eastnorth[2], double lon, double lat);
int lonlat_from_eastnorthzone(
        double out_lonlat[2],
        double e,
        double n,
        int z);

void sun_plan_point(
        int ij[2],                 // pixel coordinates in sat image frame
        double eastnorthheight[3], // utm coordinates
        double scale[2],           // satellite image resolutio
        double az, double el,      // azimuth and elevation
        double sun_height)         // sun_plan height
{
    double e = eastnorthheight[0];
    double n = eastnorthheight[1];
    double z = eastnorthheight[2];

    ij[0] = lrint( (e + (sun_height - z) * cos(el * M_PI / 180)
            * sin(az * M_PI / 180) / sin(el * M_PI / 180)) / scale[0]);

    ij[1] = lrint( - (n + (sun_height - z) * cos(el * M_PI / 180)
            * cos(az * M_PI / 180) / sin(el * M_PI / 180)) / scale[1]);
}

// get crop information for the projection of the dsm on the sun plan
void sun_plan_projection(
        int xywh[4],          // sun_plan crop information
        struct trimesh *m,    // scaled trimesh
        double scale[2],      // satellite image resolution
        double origin[3],     // relative utm origin
        double sun_height,    // sun_plan height
        double az, double el) // azimuth and elevation of the sun
{
    double xmin = INFINITY;
    double ymin = INFINITY;
    double xmax = -INFINITY;
    double ymax = -INFINITY;

    // for each vertex
    for (int i = 0; i < m->nv; i++)
    {
        // get exact utm coordinates of the vertex
        double e = m->v[3 * i + 0] + origin[0];
        double n = m->v[3 * i + 1] + origin[1];
        double z = m->v[3 * i + 2] + origin[2];

        // project this point on the sun plan
        double eastnorthheight[3] = {e, n, z};
        int ij[2] = {0, 0};
        sun_plan_point(ij, eastnorthheight, scale, az, el, sun_height);

        // check if outside the crop. If yes, enlarge the crop
        xmin = (ij[0] < xmin) ? floor(ij[0]) : xmin;
        xmax = (ij[0] > xmax) ? ceil(ij[0])  : xmax;
        ymin = (ij[1] < ymin) ? floor(ij[1]) : ymin;
        ymax = (ij[1] > ymax) ? ceil(ij[1])  : ymax;
    }

    // fill in the final crop information
    xywh[0] = lrint(xmin);
    xywh[1] = lrint(ymin);
    xywh[2] = lrint(xmax - xmin) + 1;
    xywh[3] = lrint(ymax - ymin) + 1;
}

// get direction vector of the sun rays (from sun to observer)
void sun_direction(
        double n[3],
        double az, double el)
{
    double az_rad = az * M_PI / 180;
    double el_rad = el * M_PI / 180;
    n[0] = -cos(el_rad) * sin(az_rad);
    n[0] = -cos(el_rad) * cos(az_rad);
    n[0] = -sin(el_rad);
}

// get direction vector of the camera (from camera to observer)
void camera_direction(double n[3], struct rpc *r)
{
    // give first point coordinates in satellite image
    double ijh1[3] = {500, 500, 0};

    // get this point latitude and longitude with rpc localisation
    double lonlat[2] = {0, 0};
    rpc_localization(lonlat, r, ijh1);

    // get utm coordinates from latitude and longitude
    double en[2];
    utm_from_lonlat(en, lonlat[0], lonlat[1]);

    // put utm coordinates of first point in vector director
    for (int i = 0; i < 2; i++)
        n[i] = en[i];

    // give second point coordinates in satellite image
    double ijh2[3] = {ijh1[0], ijh1[1], ijh1[2]+1};

    // get this point latitude and longitude with rpc localisation
    rpc_localization(lonlat, r, ijh2);

    // get utm coordinates from latitude and longitude
    utm_from_lonlat(en, lonlat[0], lonlat[1]);

    // create vector from the difference of the two points coordinates
    for (int i = 0; i < 2; i++)
        n[i] -= en[i];
    n[2] = ijh1[2] - ijh2[2];

    // normalise direction vector
    double norm = euclidean_norm(n, 3);
    for (int i = 0; i < 3; i++) n[i] /= norm;
}


// get crop information corresponding to the dsm projection on the sat image
void sat_im_projection(
        int  xywh[4],          // xmin, ymin, w, h of img_copy in sat image 
        struct trimesh *m,     // scaled trimesh
        struct rpc *huge_rpc,  // satellite image rpc
        double origin[3],      // relative utm origin
        int signed_zone)       // utm zone
{
    double xmin = INFINITY;
    double ymin = INFINITY;
    double xmax = -INFINITY;
    double ymax = -INFINITY;

    for (int i = 0; i < m->nv; i++)
    {
        // get utm coordinates from scaled mesh information
        double e = m->v[3 * i + 0] + origin[0];
        double n = m->v[3 * i + 1] + origin[1];
        double z = m->v[3 * i + 2] + origin[2];

        // get latitude and longitude from utm coordinates
        double lonlat[2] = {0, 0};
        lonlat_from_eastnorthzone(lonlat, e, n, signed_zone);

        // get coordinates in sat image using rpc projection
        double lonlatheight[3] = {lonlat[0], lonlat[1], z};
        double ij[2];
        rpc_projection(ij, huge_rpc, lonlatheight);

        // check if outside of current crop
        xmin = (ij[0] < xmin) ? floor(ij[0]) : xmin;
        xmax = (ij[0] > xmax) ? ceil(ij[0])  : xmax;
        ymin = (ij[1] < ymin) ? floor(ij[1]) : ymin;
        ymax = (ij[1] > ymax) ? ceil(ij[1])  : ymax;
    }

    // fill in final crop information
    xywh[0] = lrint(xmin);
    xywh[1] = lrint(ymin);
    xywh[2] = lrint(xmax) - (xmin) + 1;
    xywh[3] = lrint(ymax) - (ymin) + 1;
}

void sun_plan_pixel(
        int ij[2],                  // pixel coordinates in sun_plan
        double eastnorthheight[3],  // utm coordinates
        double scale[2],            // satellite image resolution
        double az_el[2],            // azimuth and elevation
        double sun_height,          // sun_plan height
        int xywhs[4])               // sun_plan crop information
{
        sun_plan_point(ij, eastnorthheight, scale, az_el[0],
                az_el[1], sun_height);

        for (int i = 0; i < 2; i++)
            ij[i] = ij[i] - xywhs[i];
}


void triangle_coordinates(
        double triangle_coord_meters[3][3], // utm coordinates
        int triangle_coord_pixels[3][2],    // coordinates in img_copy
        int triangle_coord_sun[3][2],       // coordinates in sun_plan
        int triangle_vertices[3],           // triangle vertices
        struct trimesh *m,                  // scaled mesh
        int xywh[4],                        // img_copy crop information
        int xywhs[4],                       // sun_plan crop information
        int signed_zone,                    // utm zone
        struct rpc *huge_rpc,               // satellite image rpc
        double origin[3],                   // relative utm origin
        double scale[2],                    // satellite image resolution
        double az_el[2],                    // azimuth and elevation
        double sun_height)                  // sun_plan height
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
            triangle_coord_meters[i][j] = m->v[3 * triangle_vertices[i] + j] + origin[j];
        double e = triangle_coord_meters[i][0]; // utm easting
        double n = triangle_coord_meters[i][1]; // utm northing
        double z = triangle_coord_meters[i][2]; // height

        // longitude and latitude from utm
        double lonlat[2] = {0, 0};
        lonlat_from_eastnorthzone(lonlat, e, n, signed_zone);
        double lonlatheight[3] = {lonlat[0], lonlat[1], z};

        // coordinates in img_copy using rpc projection
        double ij[2];
        rpc_projection(ij, huge_rpc, lonlatheight);
        for (int j = 0; j < 2; j++)
            triangle_coord_pixels[i][j] = lrint(ij[j]) - xywh[j];

        // coordinates in sun_plan using azimuth and elevation
        double eastnorthheight[3] = {e, n, z};
        sun_plan_pixel(triangle_coord_sun[i], eastnorthheight, scale, az_el,
                sun_height, xywhs);
    }
}

double determinant(double a[2], double b[2])
{
    return a[0] * b[1] - a[1] - b[0];
}

// interpolation in the case of a regular triangle
void regular_triangle_interpolation(
        int i, int j,  // pixel coordinates of the point P
        void *ee)      // triangle information
{
    struct triangle *t = ee;
    int w = t->w;

    // create the vectors PA, PB, PC, AB and AC
    int ij_a[2] = {t->coord_pixels[0][0], t->coord_pixels[0][1]};
    int ij_b[2] = {t->coord_pixels[1][0], t->coord_pixels[1][1]};
    int ij_c[2] = {t->coord_pixels[2][0], t->coord_pixels[2][1]};
    double ap[2] = {i - ij_a[0], j - ij_a[1]};
    double ab[2] = {ij_b[0] - ij_a[0], ij_b[1] - ij_a[1]};
    double ac[2] = {ij_c[0] - ij_a[0], ij_c[1] - ij_a[1]};

    // check that the determinant is different from zero (it should be if regular triangle)
    assert( determinant(ab, ac) < -0.000001 || determinant(ab, ac) > 0.000001);

    // check if highest point. If true, fill using barycentric coordinates
    if (t->img[3*(i+j*w)+2] < t->coord_meters[0][2]
                + (t->coord_meters[1][2] - t->coord_meters[0][2])
                * determinant(ap, ac)/determinant(ab, ac)
                + (t->coord_meters[2][2] - t->coord_meters[0][2])
                * determinant(ab, ap)/determinant(ab, ac))
        for (int k = 0; k < 3; k++)
            t->img[3*(i+j*w)+k]  = t->coord_meters[0][k]
                + (t->coord_meters[1][k] - t->coord_meters[0][k])
                * determinant(ap, ac)/determinant(ab, ac)
                + (t->coord_meters[2][k] - t->coord_meters[0][k])
                * determinant(ab, ap)/determinant(ab, ac);
}

// linear interpolation if the 3 points are aligned
void line_interpolation(
        int i, int j,     // coordinates of the point P
        void *ee)         // triangle information
{
    struct triangle *t = ee;
    int w = t->w;

    // create the vectors AB, AC and PA
    double ab[2] = {t->coord_pixels[1][0] - t->coord_pixels[0][0],
        t->coord_pixels[1][1] - t->coord_pixels[0][1]};
    double ac[2] = {t->coord_pixels[2][0] - t->coord_pixels[0][0],
        t->coord_pixels[2][1] - t->coord_pixels[0][1]};
    double pa[2] = {i - t->coord_pixels[0][0], j - t->coord_pixels[0][1]};

    // compute the norm of the vectors
    double ab_norm = euclidean_norm(ab, 2);
    double ac_norm = euclidean_norm(ac, 2);
    double pa_norm = euclidean_norm(pa, 2);


    if (ac_norm > 0){ // if A and C are distinct
        // check if highest point
        if ((t->coord_meters[0][2] + (t->coord_meters[2][2] -
                        t->coord_meters[0][2]) * pa_norm / ac_norm )
                > t->img[3 * (i + j * w) + 2])
            for (int k = 0; k < 3; k++)
                    t->img[3 * (i + j * w) + k] = t->coord_meters[0][k] +
                        (t->coord_meters[2][k] - t->coord_meters[0][k]) *
                        pa_norm / ac_norm;}
    else if (ab_norm > 0){ // A=C but A and B are distinct
        // check if highest point
        if ((t->coord_meters[0][2] + (t->coord_meters[1][2] -
                        t->coord_meters[1][2]) * pa_norm / ab_norm )
                > t->img[3 * (i + j * w) + 2])
            for (int k = 0; k < 3; k++)
                    t->img[3 * (i + j * w) + k] = t->coord_meters[0][k] +
                        (t->coord_meters[1][k] - t->coord_meters[0][k]) *
                        pa_norm / ab_norm;}
}

// interpolate a function in a triangle given its vertices values
void triangle_interpolation(
        int i, int j,         // pixel coordinates in image
        void *ee)             // triangle informations
{
    struct triangle *t = ee;

    // create vectors AB and AC
    double ab[2] = {t->coord_pixels[1][0] - t->coord_pixels[0][0],
        t->coord_pixels[1][1] - t->coord_pixels[0][1]};
    double ac[2] = {t->coord_pixels[2][0] - t->coord_pixels[0][0],
        t->coord_pixels[2][1] - t->coord_pixels[0][1]};

    if (determinant(ab, ac) == 0)   // if colinear (abc is a line or a point)
        line_interpolation(i, j, t);
    else                            // if non-colinear (regular triangle)
        regular_triangle_interpolation(i, j, t);
}

// fill three points if new point is highest than the previous one with
// relative utm coordinates
void fill_three_points_highest_point(
        double *img,                        // 3*w*h image
        double triangle_coord_meters[3][3], // triangle exact utm coordinates
        int triangle_coord_pixels[3][2],    // triangle pixel coordinates
        int w)                              // image width
{
    for (int i = 0; i < 3; i++)
    {
        double e = triangle_coord_meters[i][0]; // easting
        double n = triangle_coord_meters[i][1]; // northing
        double z = triangle_coord_meters[i][2]; // utm z

        int ij[2] = {triangle_coord_pixels[i][0], triangle_coord_pixels[i][1]};
        if (img[3 * (ij[1] * w + ij[0]) + 2] <= z)
        {
            img[3 * (ij[1] * w + ij[0]) + 0] = e;
            img[3 * (ij[1] * w + ij[0]) + 1] = n;
            img[3 * (ij[1] * w + ij[0]) + 2] = z;
        }
    }
}

void fill_triangle(
        double *img,                           // image to fill
        int xywh[4],                           // image size infos
        int triangle_coord_pixels[3][2],       // triangle pixel coordinates
        double triangle_coord_meters[3][3])    // triangle exact utm coordinates
{
    // create new triangle structure
    struct triangle t = {.w = xywh[2], .h = xywh[3], .img = img};

    float abc[3][2]; // vertices pixel coordinates in img
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 2; j++)
        {
            abc[i][j] = (float) triangle_coord_pixels[i][j];
            t.coord_pixels[i][j] = (float) triangle_coord_pixels[i][j];
            t.coord_meters[i][j] = triangle_coord_meters[i][j]; // utm
        }
            t.coord_meters[i][2] = triangle_coord_meters[i][2];
    }
    // fill the triangle abc using the function triangle_interpolation
    traverse_triangle(abc, triangle_interpolation, &t);
}

void initialize_img_copy_sun_plan_and_visibilities(
        double *img_copy,
        double *sun_plan,
        double *v_cam_visibility,
        double *v_sun_visibility,
        int xywh[4],
        int xywhs[4],
        int nv)
{
    // initialize sun_plan
    for (int i = 0; i < xywhs[2]*xywhs[3]; i++){
        for (int j = 0; j < 2; j++)
            sun_plan[3 * i + j] = 0;
        sun_plan[3 * i + 2] = -INFINITY;}

    // initialize img_copy
    for (int i = 0; i < xywh[2]*xywh[3]; i++){
        for (int j = 0; j < 2; j++)
            img_copy[3 * i + j] = 0;
        img_copy[3 * i + 2] = -INFINITY;}

    // initialize v_cam_visibility and v_sun_visibility
    for (int i = 0; i < nv; i++){
        v_cam_visibility[i] = false;
        v_sun_visibility[i] = false;}
}

void check_triangle_orientation(
        int t,
        double triangle_coord_meters[3][3],
        int triangle_vertices[3],
        struct trimesh *m,
        struct rpc *huge_rpc,
        double az_el[2],
        bool *v_cam_visibility,               // is the vertex seen by camera ?
        bool *v_sun_visibility,               // is the vertex seen by sun ?
        double *t_normals,                      // triangles normals
        double *t_angles)                       // triangles angles
{
        // get triangle normal and direction vectors of camera and sun
        double tri_normal[3];
        double sun_dir[3];
        double cam_dir[3];

        triangle_normal(tri_normal, triangle_coord_meters[0],
                triangle_coord_meters[1], triangle_coord_meters[2]);
        camera_direction(cam_dir, huge_rpc);
        sun_direction(sun_dir, az_el[0], az_el[1]);

        // check if triangle faces the camera
        if (scalar_product(tri_normal, cam_dir, 3) < 0)
            for (int i = 0; i < 3; i++)
                v_cam_visibility[triangle_vertices[i]] = true;

        // check if triangle faces the sun
        if (scalar_product(tri_normal, sun_dir, 3) > 0)
            for (int i = 0; i < 3; i++)
                v_sun_visibility[triangle_vertices[i]] = true;

        // compute triangle angles
        double angles[3];
        triangle_angles(angles, m, t);

        // save triangle normal and angles
        for (int i = 0; i < 3; i++){
            t_normals[3 * t + i] = tri_normal[i];
            t_angles[3 * t + i] = angles[i];}
}

// ????????????????????????????????????????
// necessite vrai triangle. Que faire pour trois points alignés distincts
static bool is_triangle(int tri_pix[3][2])
{
    return !((tri_pix[0][0] == tri_pix[1][0] && tri_pix[0][0] == tri_pix[2][0])
            || (tri_pix[0][1] == tri_pix[1][1] && tri_pix[0][1] == tri_pix[2][1]));
}


void fill_img_copy_sun_plan_with_highest_point(
        struct trimesh *m,                    // scaled mesh
        double origin[3],                     // utm origin
        int signed_zone,                      // utm zone
        struct rpc *huge_rpc,                 // rpc
        int xywh[4],                          // crop sat-im
        int xywhs[4],                         // crop sun-plan
        double scale[2],                      // image resolution
        double az_el[2],                      // azimuth and elevation
        double sun_height,                    // sun_plan height
        double *img_copy,                     // image 3*w*h
        double *sun_plan,                     // image 3*ws*hs
        bool *v_cam_visibility,               // is the vertex seen by camera ?
        bool *v_sun_visibility,               // is the vertex seen by sun ?
        double *t_normals,                      // triangles normals
        double *t_angles)                       // triangles angles
{
    for (int t = 0; t < m->nt; t++)
    {
        // get triangle coordinates in utm (with origin), pixels in the image copy and pixels in the sun plan.
        int triangle_vertices[3] = {m->t[3*t+0], m->t[3*t+1], m->t[3*t+2]};
        double triangle_coord_meters[3][3];
        int triangle_coord_pixels[3][2];
        int triangle_coord_sun[3][2];
        triangle_coordinates(triangle_coord_meters, triangle_coord_pixels,
                triangle_coord_sun, triangle_vertices, m, xywh, xywhs,
                signed_zone, huge_rpc, origin, scale, az_el, sun_height);

        // replace the projection in img_copy and sun_plan with the utm coordinates if the new point is higher than the previous one
        fill_three_points_highest_point(img_copy, triangle_coord_meters,
                triangle_coord_pixels, xywh[2]);
        fill_three_points_highest_point(sun_plan, triangle_coord_meters,
                triangle_coord_sun, xywhs[2]);

        // fill the projected triangle by interpolating the values of the vertices
        if (is_triangle(triangle_coord_pixels))
            fill_triangle(img_copy, xywh, triangle_coord_pixels, triangle_coord_meters);
        if (is_triangle(triangle_coord_sun))
            fill_triangle(sun_plan, xywhs, triangle_coord_sun, triangle_coord_meters);

        // check if triangle facing the camera and/or the sun
        check_triangle_orientation(t, triangle_coord_meters, triangle_vertices, m, huge_rpc, az_el, v_cam_visibility, v_sun_visibility, t_normals, t_angles);
    }
}

void check_vertex_visibility_satellite(
        double *out_img,                      // 3*nv image with utm info
        int v,                           // vertex index
        struct rpc *huge_rpc,                 // rpc
        int signed_zone,
        int xywh[4],                          // crop sat-im
        double scale[2],                      // image resolution
        double enz[3],                        // vertex utm coordinates
        double *img_copy,                     // image 3*w*h
        double *v_normals,                    // vertices normals
        double *v_scalars)                    // vertices normals with cam dir
{
    // get camera direction vector
    double cam_dir[3];
    camera_direction(cam_dir, huge_rpc);

    // scalar product between vertex normal and camera direction vector
    double n_vertex[3] = {v_normals[3*v], v_normals[3*v+1], v_normals[3*v+2]};
    for (int i = 0; i < 3; i++)
        v_scalars[3*v+i] = scalar_product(n_vertex, cam_dir, 3);

    // get lonlat coordinates
    double lonlat[2];
    lonlat_from_eastnorthzone(lonlat, enz[0], enz[1], signed_zone);
    double lonlatheight[3] = {lonlat[0], lonlat[1], enz[2]};

    // get coordinates in satellite image with rpc projection
    double sat_ij[2];
    rpc_projection(sat_ij, huge_rpc, lonlatheight);
    
    // get coordinates in img_copy
    int ij[2];
    for (int i = 0; i < 2; i++)
        ij[i] = lrint(sat_ij[i]) - xywh[i];

    // compute difference between vertex utm coordinates and coordinates saved in its projected pixel in img_copy
    double diff[3];
    for (int i = 0; i < 3; i++)
        diff[i] = enz[i] - img_copy[3 * (ij[1] * xywh[2] + ij[0]) + i];

    // get square scalar product between above vector and orthogonal to the direction camera vector
    double error;
    error = pow(euclidean_norm(diff,3),2) - pow(scalar_product(diff, cam_dir, 3),2);

    double sc[2] = {1.5*scale[0], 1.5*scale[1]};

    // if the vertex is close enough to the saved point, save utm coordinates
    if (img_copy[3*(ij[1]*xywh[2]+ij[0])+2] == enz[2] || pow(scalar_product(diff,cam_dir,3),2) < 4*pow(euclidean_norm(sc,2),2))
        for (int i = 0; i < 3; i++)
            out_img[3 * v + i] = enz[i];
}
    
void check_vertex_visibility_sun(
        double *out_sun,                      // 3*nv image with utm info
        int v,                                // vertex index
        int xywhs[4],                         // crop sat-im
        double scale[2],                      // image resolution
        double az_el[2],                      // sun azimuth and elevation
        double enz[3],                        // vertex utm coordinates
        double sun_height,                    // sun plan height
        double *sun_plan)                     // image 3*w*h
{
    // get sun direction vector
    double sun_dir[3];
    sun_direction(sun_dir, az_el[0], az_el[1]);

    // get coordinates in sun plan
    int ij[2];
    sun_plan_pixel(ij, enz, scale, az_el, sun_height, xywhs);

    // compute difference between vertex utm coordinates and coordinates saved in its projected pixel in img_copy
    double diff[3];
    for (int i = 0; i < 3; i++)
        diff[i] = enz[i] - sun_plan[3 * (ij[1] * xywhs[2] + ij[0]) + i];

    // get square scalar product between above vector and orthogonal to the direction camera vector
    double error;
    error = pow(euclidean_norm(diff,3), 2) - pow(scalar_product(diff, sun_dir, 3), 2);

    double sc[2] = {scale[0]+3*scale[0]/10, scale[1]+3*scale[1]/10};

    // if the vertex is close enough to the saved point, save utm coordinates
    static int counter = 0;
    static int counter2 = 0;

    static int counterA = 0;
    static int counterB = 0;
    static int counterAeB = 0;
    static int counterAoB = 0;

    bool A = (sun_plan[3*(ij[1]*xywhs[2]+ij[0])+2] == enz[2]);
    bool B = (pow(scalar_product(diff, sun_dir, 3),2) < 4*pow(euclidean_norm(sc,2),2));

    if (A){
        if (B)
            counterAeB++;
        counterA++;
    }
    if (A || B)
        counterAoB++;

    if (B)
        counterB++;
//    if ((sun_plan[3*(ij[1]*xywhs[2]+ij[0])+2] == enz[2])
//    if (sun_plan[3*(ij[1]*xywhs[2]+ij[0])+2] == enz[2]){
//             || (pow(scalar_product(diff, sun_dir, 3),2) < 4*pow(euclidean_norm(sc,2),2))){
//        if (pow(scalar_product(diff, sun_dir, 3),2) >= 4*pow(euclidean_norm(sc,2),2))
//            printf("problem diff %lf %lf %lf\n", diff[0], diff[1], diff[2]);
//        if (pow(scalar_product(diff, sun_dir, 3),2) < 4*pow(euclidean_norm(sc,2),2))
//            counter2++;
            if  (pow(scalar_product(diff, sun_dir, 3),2) < 4*pow(euclidean_norm(sc,2),2)){
//             || (pow(scalar_product(diff, sun_dir, 3),2) < 4*pow(euclidean_norm(sc,2),2))){
        counter++;
        for (int i = 0; i < 3; i++){
            if (isnan(enz[i]))
                printf("isnan v %d\n", v);
            out_sun[3 * v + i] = enz[i];
        }
    }
//    if (v > 619600)
//        printf("counter %d counter2 %d\n", counter, counter2);
    if (v > 619600)
        printf("counterA %d counterAeB %d counter B %d counterAoB %d\n", counterA, counterAeB, counterB, counterAoB);
}


void check_vertices_visibility(
        double *out_img,                      // 3*nv image with utm info
        double *out_sun,                      // 3*nv image with utm info
        struct trimesh *m,                    // scaled mesh
        double origin[3],                     // utm origin
        int signed_zone,                      // utm zone
        struct rpc *huge_rpc,                 // rpc
        int xywh[4],                          // crop sat-im
        int xywhs[4],                         // crop sun-plan
        double scale[2],                      // image resolution
        double az_el[2],                      // azimuth and elevation
        double sun_height,                    // sun_plan height
        double *img_copy,                     // image 3*w*h
        double *sun_plan,                     // image 3*ws*hs
        bool *v_cam_visibility,               // true if facing the camera
        bool *v_sun_visibility,               // true if facing the sun
        double *v_normals,                    // vertices normals
        double *v_scalars)                    // scalar product with cam dir
{
    // initialize outputs to NAN
    for (int i = 0; i < 3 * m->nv; i++){
        out_img[i] = NAN;
        out_sun[i] = NAN;}

    // initialize to vertex badly seen by satellite
    for (int i = 0; i < m->nv; i++)
        v_scalars[i] = 0;

    for (int v = 0; v < m->nv; v++)
    {
        // if the vertex faces neither sun nor camera go to next vertex
        if (!v_cam_visibility[v] && !v_sun_visibility)
           continue;

        // get utm coordinates
        double e = m->v[3 * v + 0] + origin[0];
        double n = m->v[3 * v + 1] + origin[1];
        double z = m->v[3 * v + 2] + origin[2];
        double enz[3] = {e, n, z};

        if (isnan(e) || isnan(n) || isnan(z))
            printf("enz isnan vertex v %d\n", v);

        if (v_cam_visibility[v])
            check_vertex_visibility_satellite(out_img, v, huge_rpc,
                    signed_zone, xywh, scale, enz, img_copy, v_normals, v_scalars);
        if (v_sun_visibility[v])
            check_vertex_visibility_sun(out_sun, v, xywhs, scale, az_el, enz,
                    sun_height, sun_plan);
    }
}

void add_offset_to_mesh_vertices_coordinates(
        struct trimesh *m,
        double offset[3],
        double scale[2])
{
    for (int i = 0; i < 2; i++)
        offset[i] *= scale[i];
    for (int i = 0; i < m->nv; i++)
        for (int j = 0; j < 3; j++)
            m->v[3*i+j] += offset[j];
}

int main(int c, char *v[])
{
    double ox = atof(pick_option(&c, &v, "ox", "0"));
    double oy = atof(pick_option(&c, &v, "oy", "0"));
    double oz = atof(pick_option(&c, &v, "oz", "0"));
    double xmin = atof(pick_option(&c, &v, "xmin", "0"));
    double ymin = atof(pick_option(&c, &v, "ymin", "0"));
    double offset[3] = {xmin-ox, ymin-oy, oz};
    double sun_height = atof(pick_option(&c, &v, "sh", "10000"));
    if (c < 3)
        return fprintf(stderr, "usage:\n\t"
                "%s mesh.off scale_x scale_y orig_x orig_y zone rpc"
                //0 1        2       3       4      5      6    7
                "az el utm_coord sun scalars\n", *v);
                //8 9  10        11  12
    char *filename_mesh = v[1];                     // scaled trimesh
    double scale[2] = {atof(v[2]), atof(v[3])};     // image resolution
    double origin[3] = {atof(v[4]), atof(v[5]), 0}; // utm origin
    int signed_zone = atoi(v[6]);                   // utm zone
    char *filename_rpc = v[7];                      // rpc in xml file
    double az_el[2] = {atof(v[8]), atof(v[9])};     // azimuth and elevation
    char *filename_utm_coord = v[10];               // vertices utm coordinates
    char *filename_sun = v[11];                     // vertices lightning
    char *filename_scalars = v[12];                 // vertices orientation

    // read mesh and rpc
    struct trimesh m[1];
    struct rpc huge_rpc[1];
    trimesh_read_from_off(m, filename_mesh);
    add_offset_to_mesh_vertices_coordinates(m, offset, scale);
    trimesh_fill_triangle_fans(m);
    read_rpc_file_xml(huge_rpc, filename_rpc);

    // crop informations from dsm projection on satellite image and sun
    int xywh[4];
    int xywhs[4];
    sat_im_projection(xywh, m, huge_rpc, origin, signed_zone);
    sun_plan_projection(xywhs, m, scale, origin, sun_height, az_el[0], az_el[1]);

    // images of the same size as the crop of the sun plan and satellite image
    double *sun_plan = malloc(3 * xywhs[2] * xywhs[3] * sizeof(double));
    double *img_copy = malloc(3 * xywh[2] * xywh[3] * sizeof(double));

    // vertices facing the satellite and the sun
    bool *v_cam_visibility = malloc(m->nv * sizeof(bool));
    bool *v_sun_visibility = malloc(m->nv * sizeof(bool));

    // normals and angles for triangles plus normals for angles
    double *t_normals = malloc(3 * m->nt * sizeof(double));
    double *t_angles = malloc(3 * m->nt * sizeof(double));

    // normals and scalara product of normal with camera direction vector
    double *v_normals = malloc(3 * m->nv * sizeof(double));
    double *v_scalars = malloc(3 * m->nv * sizeof(double));

    // utm coordinates for vertices if seen by satellite and/or sun
    double *out_img = malloc(3 * m->nv * sizeof(double));
    double *out_sun = malloc(3 * m->nv * sizeof(double));

    // fill img_copy and sun_plan with coordinates of highest points
    fill_img_copy_sun_plan_with_highest_point(m, origin, signed_zone, huge_rpc,
            xywh, xywhs, scale, az_el, sun_height, img_copy, sun_plan,
            v_cam_visibility, v_sun_visibility, t_normals, t_angles);

    iio_write_image_double_vec("sun_plan_get_utm.tif", sun_plan, xywhs[2], xywhs[3], 3);
    iio_write_image_double_vec("img_copy_get_utm.tif", img_copy, xywh[2], xywh[3], 3);
 
    printf("xywhs %d %d %d %d xywh %d %d %d %d", xywhs[0], xywhs[1], xywhs[2], xywhs[3], xywh[0], xywh[1], xywh[2], xywh[3]);
    // fill vertices angles
    vertices_normals_from_mesh(v_normals, m, t_angles, t_normals);

    // fill utm coordinates for vertices seen by camera and/or sun
    // get scalar product between vertex normal and camera direction vector
    check_vertices_visibility(out_img, out_sun, m, origin, signed_zone,
            huge_rpc, xywh, xywhs, scale, az_el, sun_height, img_copy,
            sun_plan, v_cam_visibility, v_sun_visibility, v_normals, v_scalars);

    // save outputs
    iio_write_image_double_vec(filename_utm_coord, out_img, m->nv, 1, 3);
    iio_write_image_double_vec(filename_sun, out_sun, m->nv, 1, 3);
    iio_write_image_double_vec(filename_scalars, v_scalars, m->nv, 1, 3);

    // free all
    free(sun_plan); free(img_copy); free(v_cam_visibility);
    free(v_sun_visibility); free(t_normals); free(t_angles); free(v_normals);
    free(v_scalars); free(out_img); free(out_sun);
}

