// create a triangular mesh from a tiff DEM (i,j,h)

#include <assert.h>   // assert
#include <math.h>     // isfinite
#include <stdio.h>    // fopen, fclose, fprintf, stdout
#include <stdlib.h>   // malloc, realloc, free, exit
#include <string.h>   // strcmp


// data structure to store a triangular mesh
struct trimesh {
	int nv, max_nv;   // number of vertices
	int nt, max_nt;   // number of triangles
	float *v;         // array of vertices
	int *t;           // array of triangles
};


// free the memory associated to the mesh
void trimesh_free_tables(struct trimesh *m)
{
	free(m->v);
	free(m->t);
}

// allocate vertex and triangle tables
static void trimesh_alloc_tables(struct trimesh *m, int nv, int nt)
{
	m->nv = 0;
	m->nt = 0;
	m->max_nv = nv;
	m->max_nt = nt;
	m->v = malloc(1000+3 * nv * sizeof*(m->v));
	m->t = malloc(1000+3 * nt * sizeof*(m->t));
}

// add a vertex (and return its index)
static int trimesh_add_vertex(struct trimesh *m, float x, float y, float z)
{
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
	return 0;
}
