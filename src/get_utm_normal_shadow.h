#ifndef GET_UTM_NORMAL_SHADOW_H
#define GET_UTM_NORMAL_SHADOW_H

#include "trimesh.h"
#include "rpc.h"

struct triangle
{
    int w, h;
    double coord_pixels[3][2];
    double coord_meters[3][3];
    double *img;
};

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
        double sun_height);        // sun_plan height

// get crop information for the projection of the dsm on the sun plan
void sun_plan_projection(
        int xywh[4],          // sun_plan crop information
        struct trimesh *m,    // scaled trimesh
        double scale[2],      // satellite image resolution
        double origin[3],     // relative utm origin
        double sun_height,    // sun_plan height
        double az, double el);// azimuth and elevation of the sun

// get direction vector of the sun rays (from sun to observer)
void sun_direction(
        double n[3],
        double az, double el);

// get direction vector of the camera (from camera to observer)
void camera_direction(double n[3], struct rpc *r);


// get crop information corresponding to the dsm projection on the sat image
void sat_im_projection(
        int  xywh[4],          // xmin, ymin, w, h of img_copy in sat image 
        struct trimesh *m,     // scaled trimesh
        struct rpc *huge_rpc,  // satellite image rpc
        double origin[3],      // relative utm origin
        int signed_zone);      // utm zone

void sun_plan_pixel(
        int ij[2],                  // pixel coordinates in sun_plan
        double eastnorthheight[3],  // utm coordinates
        double scale[2],            // satellite image resolution
        double az_el[2],            // azimuth and elevation
        double sun_height,          // sun_plan height
        int xywhs[4]);              // sun_plan crop information


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
        double sun_height);                 // sun_plan height

double determinant(double a[2], double b[2]);


// interpolation in the case of a regular triangle
void regular_triangle_interpolation(
        int i, int j,  // pixel coordinates of the point P
        void *ee);     // triangle information

// linear interpolation if the 3 points are aligned
void line_interpolation(
        int i, int j,     // coordinates of the point P
        void *ee);        // triangle information

// interpolate a function in a triangle given its vertices values
void triangle_interpolation(
        int i, int j,         // pixel coordinates in image
        void *ee);            // triangle informations 

// fill three points if new point is highest than the previous one with 
// relative utm coordinates
void fill_three_points_highest_point(
        double *img,                        // 3*w*h image
        double triangle_coord_meters[3][3], // triangle exact utm coordinates 
        int triangle_coord_pixels[3][2],    // triangle pixel coordinates
        int w);                             // image width

void fill_triangle(
        double *img,                           // image to fill
        int xywh[4],                           // image size infos
        int triangle_coord_pixels[3][2],       // triangle pixel coordinates
        double triangle_coord_meters[3][3]);   // triangle exact utm coordinates

void initialize_img_copy_sun_plan_and_visibilities(
        double *img_copy,
        double *sun_plan,
        double *v_cam_visibility,
        double *v_sun_visibility,
        int xywh[4], 
        int xywhs[4], 
        int nv);

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
        double *t_angles);                      // triangles angles

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
        double *t_angles);                      // triangles angles

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
        double *v_scalars);                   // vertices normals with cam dir
    
void check_vertex_visibility_sun(
        double *out_sun,                      // 3*nv image with utm info
        int v,                                // vertex index
        int xywhs[4],                         // crop sat-im
        double scale[2],                      // image resolution
        double az_el[2],                      // sun azimuth and elevation
        double enz[3],                        // vertex utm coordinates
        double sun_height,                    // sun plan height
        double *sun_plan);                    // image 3*w*h
    

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
        double *v_scalars);                    // scalar product with cam dir 

void add_offset_to_mesh_vertices_coordinates(
        struct trimesh *m, 
        double offset[3],
        double scale[2]);

#endif
