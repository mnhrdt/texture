#ifndef TRIMESH_H
#define TRIMESH_H

// data structure and functions for triangular mesh I/O
// #defines
#define TRIMESH_MORE_STUFF


// "struct trimesh" : a data structure to store a triangular mesh
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


// free the memory associated to the mesh
void trimesh_free_tables(struct trimesh *m);

// function to create a mesh from a digital elevation map
void trimesh_create_from_dem(
        struct trimesh *m,
        float *x,
        int w, int h);

// function to create a mesh from a digital elevation map with offset
void trimesh_create_from_dem_with_scale(
        struct trimesh *m,
        float *x,
        int w,
        int h,
        double scale[3]);

// function to save a triangulated surface to a ply file
void trimesh_write_to_ply(
        char *fname,
        struct trimesh *m);

// function to save a cloud point to a coloured ply file
void cloud_write_to_coloured_ply(
        char *fname,
        struct trimesh *m,
        double *c,
        double origin[3]);

void trimesh_write_to_coloured_ply(char *fname,
	       	struct trimesh *m, double *c, double origin[3]);
void trimesh_write_to_off(char *fname, struct trimesh *m);
void trimesh_write_to_coloured_off(char *fname,
		struct trimesh *m, double *c, double origin[3]);

void trimesh_read_from_ply(struct trimesh *m, char *fname);

void trimesh_read_from_off(struct trimesh *m, char *fname);

#ifdef TRIMESH_MORE_STUFF
int uniq(void *X, int n, int s, int (*c)(const void *, const void *));
void trimesh_fill_edges(struct trimesh *m);
void trimesh_fill_triangle_fans(struct trimesh *m);
int trimesh_get_triangle_fan(int *out, struct trimesh *m, int v);
void trimesh_dump_edges(char *filename, struct trimesh *m);
#endif//TRIMESH_MORE_STUFF


#ifdef TRIMESH_CALCULUS
#error "clarify double edges stuff"
void trimesh_gradient(struct trimesh *m, float *y, float *x);
void trimesh_centering(struct trimesh *m, float *y, float *x);
void trimesh_recentering( struct trimesh *m, float *y, float *x);
void trimesh_divergence(struct trimesh *m, float *y, float *x);
void trimesh_laplacian(struct trimesh *m, float *y, float *x);
#endif//TRIMESH_CALCULUS
#endif
