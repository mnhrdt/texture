// data structure and functions for triangular mesh I/O



#include <assert.h>   // assert
#include <math.h>     // isfinite
#include <stdbool.h>  // bool
#include <stdio.h>    // fopen, fclose, fprintf, stdout
#include <stdlib.h>   // malloc, realloc, free, exit
#include <string.h>   // strcmp


// data structure to store a triangular mesh
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

	// 2.1. graph of edges (a directed graph with double edges)
	int ne;       // number of edges
	int *e;       // list of edges (pairs of vertex indices)

	// 2.2. access from any point
	// sides of edges (ne pairs of triangle indices)
	// sides of triangles (nt triplets of edge indices)

	// 2.3. edge fans (two tables of length nv)
	// first edge around each vertex
	// next edge in the triangle list

	// 2.4. triangle fans
	// first triangle around each vertex (two tables of length nv)
	// next triangle in the triangle list

        // 3. color data that has to be filled during colorisation.
        float *vc; // vertices (RGB)
#endif//TRIMESH_MORE_STUFF
};


// free the memory associated to the mesh
void trimesh_free_tables(struct trimesh *m)
{
	free(m->v);
	free(m->t);
}

// allocate vertex and triangle tables
// (this function is always called to construcit a mesh)
static void trimesh_alloc_tables(struct trimesh *m, int nv, int nt)
{
	m->nv = 0;
	m->nt = 0;
	m->max_nv = nv;
	m->max_nt = nt;
	m->v = malloc(1000+3 * nv * sizeof*(m->v));
	m->t = malloc(1000+3 * nt * sizeof*(m->t));

#ifdef TRIMESH_MORE_STUFF
	m->ne = 0;
	m->e = NULL;
        m->vc = NULL;
#endif//TRIMESH_MORE_STUFF
}

// add a vertex (and return its index)
static int trimesh_add_vertex(struct trimesh *m, float x, float y, float z)
{
	//fprintf(stderr, "trimesh_add_vertex_%d : %g %g %g\n", m->nv, x,y,z);
	// TODO: instead of aborting, realloc more memory as necessary
	assert(m->nv + 1 <= m->max_nv);

	m->v[3*m->nv + 0] = x;
	m->v[3*m->nv + 1] = y;
	m->v[3*m->nv + 2] = z;
	return m->nv++;
}

// add a triangle (and return its index)
static int trimesh_add_triangle(struct trimesh *m, int a, int b, int c)
{
	//fprintf(stderr, "trimesh_add_triangle_%d : %d %d %d\n", m->nt, a,b,c);
	// TODO: instead of aborting, realloc more memory as necessary
	assert(m->nt + 1 <= m->max_nt);

	assert(a >= 0); assert(a < m->nv);
	assert(b >= 0); assert(b < m->nv);
	assert(c >= 0); assert(c < m->nv);

	m->t[3*m->nt + 0] = a;
	m->t[3*m->nt + 1] = b;
	m->t[3*m->nt + 2] = c;
	return m->nt++;
}

// function to create a mesh from a digital elevation map
void trimesh_create_from_dem(struct trimesh *m, float *x, int w, int h)
{
	// initialize the mesh
	trimesh_alloc_tables(m, w*h, 2*(w-1)*(h-1));

	// build a table of vertex correspondences
	int *t = malloc(w * h * sizeof*t);
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	if (isfinite(x[j*w + i]))
		t[j*w + i] = trimesh_add_vertex(m, i, j, x[j*w + i]);
	else
		t[j*w + i] = -1;

	// add the triangles all of whose whose vertices are good
	for (int j = 0; j < h-1; j++)
	for (int i = 0; i < w-1; i++)
	{
		int a = (j+0)*w + (i+0);
		int b = (j+0)*w + (i+1);
		int c = (j+1)*w + (i+0);
		int d = (j+1)*w + (i+1);

		// TODO: criteria to choose the appropriate diagonal
		if (t[a] >= 0 && t[b] >= 0 && t[c] >= 0)
			trimesh_add_triangle(m, t[a], t[b], t[c]);
		if (t[c] >= 0 && t[b] >= 0 && t[d] >= 0)
			trimesh_add_triangle(m, t[c], t[b], t[d]);
	}

	// cleanup
	free(t);
}

// function to save a triangulated surface to a ply file
void trimesh_write_to_ply(char *fname, struct trimesh *m)
{
	// dump the ply file (with dsm-inherited connectivity)
	FILE *f = strcmp(fname, "-") ? fopen(fname, "w") : stdout;
	if (!f)
		exit(fprintf(stderr, "ERROR: cannot open file (%s)\n", fname));

	// print header
	fprintf(f, "ply\n");
	fprintf(f, "format ascii 1.0\n");
	fprintf(f, "comment bare triangulated surface\n");
	fprintf(f, "element vertex %d\n", m->nv);
	fprintf(f, "property float x\n");
	fprintf(f, "property float y\n");
	fprintf(f, "property float z\n");
	fprintf(f, "element face %d\n", m->nt);
	fprintf(f, "property list uchar int vertex_indices\n");
	fprintf(f, "end_header\n");

	// print points
	for (int i = 0; i < m->nv; i++)
                fprintf(f, "%.16lf %.16lf %.16lf\n",
				m->v[3*i+0], m->v[3*i+1], m->v[3*i+2]);

	// print triangles
	for (int i = 0; i < m->nt; i++)
                fprintf(f, "3 %d %d %d\n",
				m->t[3*i+0], m->t[3*i+1], m->t[3*i+2]);

	// cleanup
	fclose(f);
}

void trimesh_read_from_ply(struct trimesh *m, char *fname)
{
	FILE *f = strcmp(fname, "-") ? fopen(fname, "r") : stdin;
	if (!f)
		exit(fprintf(stderr, "ERROR: cannot open file (%s)\n", fname));

	int n_vertices = -1;
	int n_triangles = -1;

	// process header lines
	char buf[FILENAME_MAX] = {0};
	while (fgets(buf, FILENAME_MAX, f)) { // read a line into "buf"
		int tmp;
		if (1 == sscanf(buf, "element vertex %d\n", &tmp))
			n_vertices = tmp;
		if (1 == sscanf(buf, "element face %d\n", &tmp))
			n_triangles = tmp;
		if (0 == strcmp(buf, "end_header\n"))
			break;
	}
	//fprintf(stderr, "n_vertices = %d\n", n_vertices);
	//fprintf(stderr, "n_triangles = %d\n", n_triangles);

	// create mesh structure
	trimesh_alloc_tables(m, n_vertices, n_triangles);

	// read vertices
	while (m->nv < n_vertices && fgets(buf, FILENAME_MAX, f))
	{
		double x[3];
		if (3 != sscanf(buf, "%lf %lf %lf\n", x, x+1, x+2))
			exit(fprintf(stderr, "ERROR: vfail_1 \"%s\"\n", buf));
		trimesh_add_vertex(m, x[0], x[1], x[2]);
	}
	if (n_vertices != m->nv)
		exit(fprintf(stderr, "ERROR: nv %d %d\n", n_vertices, m->nv));

	// read triangles
	while (m->nt < n_triangles && fgets(buf, FILENAME_MAX, f))
	{
		int x[3];
		if (3 != sscanf(buf, "3 %d %d %d\n", x, x+1, x+2))
			exit(fprintf(stderr, "ERROR: tfail_1 \"%s\"\n", buf));
		trimesh_add_triangle(m, x[0], x[1], x[2]);
	}
	if (n_triangles != m->nt)
		exit(fprintf(stderr, "ERROR: nt %d %d\n", n_triangles, m->nt));

	// cleanup and exit
	fclose(f);
}

#ifdef TRIMESH_MORE_STUFF

#define BAD_MIN(a,b) (a)<(b)?(a):(b);
#define BAD_MAX(a,b) (a)>(b)?(a):(b);
void trimesh_fill_edges(struct trimesh *m)
{
	// allocate space for the maximum possible number of edges
	int *e = malloc(2 * 3 * m->nt * sizeof*e);

	// add all edges in a canonical orientation, possibly repeated
	int ne = 0;
	for (int i = 0; i < m->nt; i++)
	for (int k = 0; k < 3; k++)
	{
		e[2*ne+0] = BAD_MIN(m->t[3*i+k], m->t[3*i+(k%3)]);
		e[2*ne+1] = BAD_MAX(m->t[3*i+k], m->t[3*i+(k%3)]);
		ne += 1;
	}

	// sort the list of edges
	qsort(e, ne, 2*sizeof*e, compare_int_pair);

	// remove repeated edges
	int *p = e;
	int *q = e;
	while (q < e + 2*ne)
	{
		if (p[0] != q[0] || p[1] != q[1])
		{
			p[0] = q[0];
			p[1] = q[1];
			p += 2;
		}
		q += 2;
	}

	// duplicate the list of edges
	for (int i = 0; i < ne; i++)
	{
		e[2*(ne+i) + 0] = e[2*i + 1];
		e[2*(ne+i) + 1] = e[2*i + 0];
	}

	// update struct, and finish
	m->e = e;
	m->ne = 2 * ne;
}
#endif//TRIMESH_MORE_STUFF


#ifdef TRIMESH_CALCULUS
void trimesh_gradient(
		struct trimesh *m,  // base mesh
		float *y,           // output gradient (array of length m->ne/2)
		float *x            // input function (array of length m->nv)
		)
{
	for (int i = 0; i < m->ne/2; i++)
		y[i] = x[ m->e[2*i+1] ]  -  x[ m->e[2*i+1] ];
}

void trimesh_centering(
		struct trimesh *m,  // base mesh
		float *y,           // output field (array of length m->ne/2)
		float *x            // input function (array of length m->nv)
		)
{
	for (int i = 0; i < m->ne/2; i++)
		y[i] = 0.5 * (  x[ m->e[2*i+1] ]  +  x[ m->e[2*i+1] ]  );
}

void trimesh_divergence(
		struct trimesh *m,  // base mesh
		float *y,           // output divergence (array of length m->nv)
		float *x            // input field (array of length m->ne/2)
		)
{
	for (int i = 0; i < m->nv; i++)
		y[i] = 0;

	for (int i = 0; i < m->ne/2; i++)
	{
		y[ m->e[2*i+1] ] += x[i];
		y[ m->e[2*i+0] ] -= x[i];
	}
}

void trimesh_laplacian(
		struct trimesh *m,  // base mesh
		float *y,           // output laplacian (array of length m->nv)
		float *x            // input function (array of length m->nv)
		)
{
	float *d = malloc(m->ne * sizeof*d);
	trimesh_gradient(m, d, x);
	trimesh_divergence(m, y, d);
	free(d);
}
#endif//TRIMESH_CALCULUS

#ifdef TRIMESH_DEMO_MAIN
#include "iio.h"
int main(int c, char *v[])
{
	if (c != 2 && c != 3)
		return fprintf(stderr, "usage:\n\t%s in.tif [out.ply]\n", *v);
	//                                         0 1       2
	char *filename_in  = v[1];
	char *filename_out = c > 2 ? v[2] : "-";

	// read input DSM
	int w, h;
	float *x = iio_read_image_float(filename_in, &w, &h);

	// create triangulation
	struct trimesh m[1];
	trimesh_create_from_dem(m, x, w, h);

	// save triangulation into ply flie
	trimesh_write_to_ply(filename_out, m);


	// cleanup and exit
	trimesh_free_tables(m);
	free(x);

	// do a silly consistency loop
	struct trimesh n[1];
	trimesh_read_from_ply(n, filename_out);
	trimesh_write_to_ply("/tmp/plytmp.ply", n);

	// cleanup and exit
	trimesh_free_tables(n);
	//free(x);
	return 0;
}
#endif//TRIMESH_DEMO_MAIN
