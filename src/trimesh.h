#ifndef TRIMESH_H
#define TRIMESH_H

// data structure and functions for triangular mesh I/O                     {{{1
// #defines                                                                 {{{1
#define TRIMESH_MORE_STUFF


// "struct trimesh" : a data structure to store a triangular mesh           {{{1
struct trimesh {
	// 1. essential data
	int nv;       // number of vertices
	int nt;       // number of triangles
	int max_nv;   // size of vertex array
	int max_nt;   // size of triangle array
	float *v;     // vertices (3D points)
	int *t;       // triangles (triplets of vertex indices)

#ifdef TRIMESH_MORE_STUFF
	// 2. accessory data, that can be computed from the essential data

	// 2.1. graph of edges (a directed graph, without double edges)
	int ne;       // number of edges
	int *e;       // list of edges (pairs of vertex indices)

	// 2.2. access from any point
	// sides of edges (ne pairs of triangle indices)
	// sides of triangles (nt triplets of edge indices)

	// 2.3. edge fans (two tables of length nv)
	// first edge around each vertex
	// next edge in the edge list

	// 2.4. triangle fans
	int *tfirst;  // first triangle around each vertex (length nv)
	int *tnext_a; // next triangle in the triangle list (length nt)
	int *tnext_b; // next triangle in the triangle list (length nt)
	int *tnext_c; // next triangle in the triangle list (length nt)

#endif//TRIMESH_MORE_STUFF
};


// free the memory associated to the mesh                                   {{{1
void trimesh_free_tables(struct trimesh *m);

// function to create a mesh from a digital elevation map                   {{{1
void trimesh_create_from_dem(
        struct trimesh *m, 
        float *x, 
        int w, int h);

// function to create a mesh from a digital elevation map with offset       {{{1
void trimesh_create_from_dem_with_scale(
        struct trimesh *m, 
        float *x, 
        int w, 
        int h,
        double scale[3]);

// function to save a triangulated surface to a ply file                    {{{1
void trimesh_write_to_ply(
        char *fname, 
        struct trimesh *m);

// function to save a cloud point to a coloured ply file           {{{1
void cloud_write_to_coloured_ply(
        char *fname, 
        struct trimesh *m, 
        double *c, 
        double origin[3]);

// function to save a cloud point to a coloured ply file           {{{1
void trimesh_write_to_coloured_ply(
        char *fname, 
        struct trimesh *m, 
        double *c, 
        double origin[3]);

// function to save a triangulated surface to an off file                    {{{1
void trimesh_write_to_off(
        char *fname, 
        struct trimesh *m);

// function to save a triangulated surface to an off file                    {{{1
void trimesh_write_to_coloured_off(
        char *fname,          // where to save the mesh
        struct trimesh *m,    // scaled mesh with relative utm coords
        double *c,            // vertices rgb (int between 0 and 255)
        double origin[3]);    // utm origin

// function to read a triangulated surface from an ply file                  {{{1
void trimesh_read_from_ply(
        struct trimesh *m, 
        char *fname);

// function to read a triangulated surface from a off file                  {{{1
void trimesh_read_from_off(
        struct trimesh *m, 
        char *fname);
#ifdef TRIMESH_MORE_STUFF

// function to compute the array of edges of a mesh                         {{{1
#define BAD_MIN(a,b) (a)<(b)?(a):(b);
#define BAD_MAX(a,b) (a)>(b)?(a):(b);
int uniq(
        void *X, 
        int n, 
        int s, 
        int (*c)(const void *, const void *));

void trimesh_fill_edges(struct trimesh *m);

void trimesh_fill_triangle_fans(struct trimesh *m);

int trimesh_get_triangle_fan(int *out, struct trimesh *m, int v);

#endif//TRIMESH_MORE_STUFF


#ifdef TRIMESH_CALCULUS
#error "clarify double edges stuff"
void trimesh_gradient(                                                   // {{{1
		struct trimesh *m,  // base mesh
		float *y,           // output gradient (array of length m->ne/2)
		float *x            // input function (array of length m->nv)
		);

void trimesh_centering(                                                  // {{{1
		struct trimesh *m,  // base mesh
		float *y,           // output field (array of length m->ne/2)
		float *x            // input function (array of length m->nv)
		);

void trimesh_recentering(                                                // {{{1
		struct trimesh *m,  // base mesh
		float *y,           // output function (array of length m->nv)
		float *x            // input field (array of length m->ne/2)
		);

void trimesh_divergence(                                                 // {{{1
		struct trimesh *m,  // base mesh
		float *y,           // output divergence (array of length m->nv)
		float *x            // input field (array of length m->ne/2)
		);

void trimesh_laplacian(                                                  // {{{1
		struct trimesh *m,  // base mesh
		float *y,           // output laplacian (array of length m->nv)
		float *x            // input function (array of length m->nv)
		);
#endif//TRIMESH_CALCULUS
#endif

// vim:set foldmethod=marker:
