// data structure and functions for triangular mesh I/O                     {{{1

// #includes                                                                {{{1

#include <assert.h>   // assert
#include <math.h>     // isfinite
#include <stdbool.h>  // bool
#include <stdio.h>    // fopen, fclose, fprintf, stdout
#include <stdlib.h>   // malloc, realloc, free, exit
#include <string.h>   // strcmp

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

	// 2.1. graph of edges (a directed graph with double edges)
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
void trimesh_free_tables(struct trimesh *m)
{
	free(m->v);
	free(m->t);
#ifdef TRIMESH_MORE_STUFF
	void *t[5] = {m->e, m->tfirst, m->tnext_a, m->tnext_b, m->tnext_c};
	for (int i = 0; i < 5; i++)
		if (t[i])
			free(t[i]);
#endif//TRIMESH_MORE_STUFF
}

// allocate vertex and triangle tables                                      {{{1
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
	m->tfirst = m->tnext_a = m->tnext_b = m->tnext_c = NULL;
#endif//TRIMESH_MORE_STUFF
}

// add a vertex (and return its index)                                      {{{1
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

// add a triangle (and return its index)                                    {{{1
static void swapint(int *a, int *b)
{
	int t = *a;
	*a = *b;
	*b = t;
}
static void sort_three_ints(int *y, int *x)
{
	y[0]=x[0];
	y[1]=x[1];
	y[2]=x[2];
	if (y[0] > y[1]) swapint(y+0, y+1);
	if (y[1] > y[2]) swapint(y+1, y+2);
	if (y[0] > y[1]) swapint(y+0, y+1);
}
static int trimesh_add_triangle(struct trimesh *m, int a, int b, int c)
{
	//fprintf(stderr, "trimesh_add_triangle_%d : %d %d %d\n", m->nt, a,b,c);
	// TODO: instead of aborting, realloc more memory as necessary
	assert(m->nt + 1 <= m->max_nt);

	assert(a >= 0); assert(a < m->nv);
	assert(b >= 0); assert(b < m->nv);
	assert(c >= 0); assert(c < m->nv);

	// XXX : the following three lines must remain commented,
	// otherwise, the orientation of the mesh is lost
	//if (a > b) swapint(&a, &b);
	//if (b > c) swapint(&b, &c);
	//if (a > b) swapint(&a, &b);

	m->t[3*m->nt + 0] = a;
	m->t[3*m->nt + 1] = b;
	m->t[3*m->nt + 2] = c;
	return m->nt++;
}

// function to create a mesh from a digital elevation map                   {{{1
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

// function to create a mesh from a digital elevation map with offset       {{{1
void trimesh_create_from_dem_with_scale(struct trimesh *m, 
        float *x, 
        int w, 
        int h,
        double scale[3])
{
       // initialize the mesh
       trimesh_alloc_tables(m, w*h, 2*(w-1)*(h-1));

       // build a table of vertex correspondences
       int *t = malloc(w * h * sizeof*t);
       for (int j = 0; j < h; j++)
       for (int i = 0; i < w; i++)
       if (isfinite(x[(j*w+i)%(w*h)]))
               t[(j*w+i)%(w*h)] = trimesh_add_vertex(m, 
                       scale[0] * i,
                       scale[1] * j,
                       scale[2] * x[(j*w+i)%(w*h)]);
       else
               t[(j*w+i)%(w*h)] = -1;

       // add the triangles all of those whose vertices are good
       for (int j = 0; j < h-1; j++)
       for (int i = 0; i < w-1; i++)
       {
               int a = ((j+0)*w + (i+0))%(w*h);
               int b = ((j+0)*w + (i+1))%(w*h);
               int c = ((j+1)*w + (i+0))%(w*h);
               int d = ((j+1)*w + (i+1))%(w*h);

                if (fabs(x[b]-x[c]) < fabs(x[a]-x[d]))
                {
               if (t[a] >= 0 && t[b] >= 0 && t[c] >= 0)
                       trimesh_add_triangle(m, t[a], t[b], t[c]);
               if (t[c] >= 0 && t[b] >= 0 && t[d] >= 0)
                       trimesh_add_triangle(m, t[c], t[b], t[d]);
                }
                else
                {
               if (t[a] >= 0 && t[b] >= 0 && t[d] >= 0)
                       trimesh_add_triangle(m, t[a], t[b], t[d]);
               if (t[a] >= 0 && t[c] >= 0 && t[d] >= 0)
                       trimesh_add_triangle(m, t[a], t[d], t[c]);
                }

       }

}
// function to create a mesh from a digital elevation map with offset       {{{1
void trimesh_create_from_dem_with_offset(struct trimesh *m, 
        float *x, 
        int w, 
        int h,
        double ox,
        double oy,
        double oz)
{
       // initialize the mesh
       trimesh_alloc_tables(m, w*h, 2*(w-1)*(h-1));

       // build a table of vertex correspondences
       int *t = malloc(w * h * sizeof*t);
       for (int j = oy; j < h+oy; j++)
       for (int i = ox; i < w+ox; i++)
       if (isfinite(x[(j*w+i)%(w*h)]))
               t[(j*w+i)%(w*h)] = trimesh_add_vertex(m, i, j, x[(j*w+i)%(w*h)]);
       else
               t[(j*w+i)%(w*h)] = -1;

       // add the triangles all of those whose vertices are good
       for (int j = oy; j < h-1+oy; j++)
       for (int i = ox; i < w-1+ox; i++)
       {
               int a = ((j+0)*w + (i+0))%(w*h);
               int b = ((j+0)*w + (i+1))%(w*h);
               int c = ((j+1)*w + (i+0))%(w*h);
               int d = ((j+1)*w + (i+1))%(w*h);

                if (fabs(x[b]-x[c]) < fabs(x[a]-x[d]))
                {
               if (t[a] >= 0 && t[b] >= 0 && t[c] >= 0)
                       trimesh_add_triangle(m, t[a], t[b], t[c]);
               if (t[c] >= 0 && t[b] >= 0 && t[d] >= 0)
                       trimesh_add_triangle(m, t[c], t[b], t[d]);
                }
                else
                {
               if (t[a] >= 0 && t[b] >= 0 && t[d] >= 0)
                       trimesh_add_triangle(m, t[a], t[b], t[d]);
               if (t[a] >= 0 && t[c] >= 0 && t[d] >= 0)
                       trimesh_add_triangle(m, t[a], t[d], t[c]);
                }

       }

       // cleanup
       free(t);
}
// function to save a triangulated surface to a ply file                    {{{1
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

// function to save a triangulated surface to an off file with offset        {{{1
void trimesh_write_to_off_with_offset(char *fname, struct trimesh *m, 
        double ox, double oy, double oz)
{
	// dump the off file (with dsm-inherited connectivity)
	FILE *f = strcmp(fname, "-") ? fopen(fname, "w") : stdout;
	if (!f)
		exit(fprintf(stderr, "ERROR: cannot open file (%s)\n", fname));

	// print header
	fprintf(f, "OFF\n");
	fprintf(f, "%d %d 0\n", m->nv, m->nt);

	// print points
	for (int i = 0; i < m->nv; i++)
                fprintf(f, "%.16lf %.16lf %.16lf\n",
				m->v[3*i+0] + ox, 
                                m->v[3*i+1] + oy, 
                                m->v[3*i+2] + oz);

	// print triangles
	for (int i = 0; i < m->nt; i++)
                fprintf(f, "3 %d %d %d\n",
				m->t[3*i+0], m->t[3*i+1], m->t[3*i+2]);

	// cleanup
	fclose(f);
}
// function to save a triangulated surface to an off file                    {{{1
void trimesh_write_to_off(char *fname, struct trimesh *m)
{
	// dump the off file (with dsm-inherited connectivity)
	FILE *f = strcmp(fname, "-") ? fopen(fname, "w") : stdout;
	if (!f)
		exit(fprintf(stderr, "ERROR: cannot open file (%s)\n", fname));

	// print header
	fprintf(f, "OFF\n");
	fprintf(f, "%d %d 0\n", m->nv, m->nt);

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
// function to save a triangulated surface to an off file                    {{{1
// with dimensions in meters
void trimesh_write_to_scaled_off(char *fname, struct trimesh *m)
{
	// dump the off file (with dsm-inherited connectivity)
	FILE *f = strcmp(fname, "-") ? fopen(fname, "w") : stdout;
	if (!f)
		exit(fprintf(stderr, "ERROR: cannot open file (%s)\n", fname));

	// print header
	fprintf(f, "OFF\n");
	fprintf(f, "%d %d 0\n", m->nv, m->nt);

	// print points
	for (int i = 0; i < m->nv; i++)
                fprintf(f, "%.16lf %.16lf %.16lf\n",
				0.3*m->v[3*i+0], 0.3*m->v[3*i+1], m->v[3*i+2]);

	// print triangles
	for (int i = 0; i < m->nt; i++)
                fprintf(f, "3 %d %d %d\n",
				m->t[3*i+0], m->t[3*i+1], m->t[3*i+2]);

	// cleanup
	fclose(f);
}
// function to save a triangulated surface to a coloured ply file           {{{1
void trimesh_write_to_coloured_ply(char *fname, struct trimesh *m, double *c, double t)
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
        fprintf(f, "property uchar red\n");
        fprintf(f, "property uchar green\n");
        fprintf(f, "property uchar blue\n");
	fprintf(f, "element face %d\n", m->nt);
	fprintf(f, "property list uchar int vertex_indices\n");
	fprintf(f, "end_header\n");
        

	// print points
        for (int i = 0; i < m->nv; i++)
        {
            if (isnan(c[3*i+0]))
                fprintf(f, "%.16lf %.16lf %.16lf 0 0 255\n",
                        0.3*m->v[3*i+0], -0.3*m->v[3*i+1],
                        m->v[3*i+2]);
            else
                fprintf(f, "%.16lf %.16lf %.16lf %d %d %d\n",
                        0.3*m->v[3*i+0], -0.3*m->v[3*i+1],
                        m->v[3*i+2],
                        (int) fmax(fmin(c[3*i+0],255),0),
                        (int) fmax(fmin(c[3*i+1],255),0),
                        (int) fmax(fmin(c[3*i+2],255),0));
        }


	// print triangles
	for (int i = 0; i < m->nt; i++)
                fprintf(f, "3 %d %d %d\n",
				m->t[3*i+0], m->t[3*i+1], m->t[3*i+2]);

	// cleanup
	fclose(f);
}

// function to save a triangulated surface to a coloured ply file           {{{1
void trimesh_write_to_coloured_ply2(char *fname, struct trimesh *m, 
        double *c, double origin[3])
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
        fprintf(f, "property uchar red\n");
        fprintf(f, "property uchar green\n");
        fprintf(f, "property uchar blue\n");
	fprintf(f, "element face %d\n", m->nt);
	fprintf(f, "property list uchar int vertex_indices\n");
	fprintf(f, "end_header\n");
        

	// print points
        for (int i = 0; i < m->nv; i++)
        {
            if (isnan(c[3*i+0]))
                fprintf(f, "%.16lf %.16lf %.16lf 0 0 255\n",
                        m->v[3*i+0] + origin[0], 
                        m->v[3*i+1] + origin[1],
                        m->v[3*i+2] + origin[2]);
            else
                fprintf(f, "%.16lf %.16lf %.16lf %d %d %d\n",
                        m->v[3*i+0] + origin[0], 
                        m->v[3*i+1] + origin[1], 
                        m->v[3*i+2] + origin[2], 
                        (int) fmax(fmin(c[3*i+0],255),0),
                        (int) fmax(fmin(c[3*i+1],255),0),
                        (int) fmax(fmin(c[3*i+2],255),0));
        }


	// print triangles
	for (int i = 0; i < m->nt; i++)
                fprintf(f, "3 %d %d %d\n",
				m->t[3*i+0], m->t[3*i+1], m->t[3*i+2]);

	// cleanup
	fclose(f);
}

// function to read a triangulated surface from a off file                  {{{1
void trimesh_read_from_off(struct trimesh *m, char *fname)
{
	FILE *f = strcmp(fname, "-") ? fopen(fname, "r") : stdin;
	if (!f)
		exit(fprintf(stderr, "ERROR: cannot open file (%s)\n", fname));

	int n_vertices = -1;
	int n_triangles = -1;

	// process header lines
	char buf[FILENAME_MAX] = {0};

        if (2 != fscanf(f, "OFF\n%d %d 0\n", &n_vertices, &n_triangles))
            exit(fprintf(stderr, "ERROR: cannot read number of vertices and triangles\n"));
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

// function to read a triangulated surface from a off file                  {{{1
void trimesh_read_from_scaled_off(struct trimesh *m, char *fname)
{
	FILE *f = strcmp(fname, "-") ? fopen(fname, "r") : stdin;
	if (!f)
		exit(fprintf(stderr, "ERROR: cannot open file (%s)\n", fname));

	int n_vertices = -1;
	int n_triangles = -1;

	// process header lines
	char buf[FILENAME_MAX] = {0};

        if (2 != fscanf(f, "OFF\n%d %d 0\n", &n_vertices, &n_triangles))
            exit(fprintf(stderr, "ERROR: cannot read number of vertices and triangles\n"));
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
		trimesh_add_vertex(m, x[0]/0.3, x[1]/0.3, x[2]);
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

// function to read a triangulated surface from an ply file                  {{{1
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

// function to compute the array of edges of a mesh                         {{{1
#define BAD_MIN(a,b) (a)<(b)?(a):(b);
#define BAD_MAX(a,b) (a)>(b)?(a):(b);
static int compare_int_pair(const void *aa, const void *bb)
{
	const int *a = (const int *)aa;
	const int *b = (const int *)bb;
	int x = (a[0] < b[0]) - (a[0] > b[0]);
	int y = (a[1] < b[1]) - (a[1] > b[1]);
	return x ? x : y;
}
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

static void add_triangle_to_list(int *first, int *next, int t, int v, char *s)
{
	// update a linked list (with -1 to indicate the end of each chains)
	if (first[v] < 0)
		first[v] = t;
	else
		next[t] = first[v];
	first[v] = t;
}

// function to compute the triangle fans of a mesh                         {{{1
static void trimesh_fill_triangle_fans(struct trimesh *m)
{
	// allocate space
	int *f = malloc(m->nv * sizeof*f);
	int *a = malloc(m->nt * sizeof*a);
	int *b = malloc(m->nt * sizeof*b);
	int *c = malloc(m->nt * sizeof*c);

	// initialize each linked list to empty
	for (int i = 0; i < m->nv; i++) f[i] = -1;
	for (int i = 0; i < m->nt; i++) a[i] = -1;
	for (int i = 0; i < m->nt; i++) b[i] = -1;
	for (int i = 0; i < m->nt; i++) c[i] = -1;

	// add each triangle to its list
	for (int i=0; i < m->nt; i++)
	{
		int abc[3];
		sort_three_ints(abc, m->t + 3*i);
		add_triangle_to_list(f, a, i, abc[0], "a");
		add_triangle_to_list(f, b, i, abc[1], "b");
		add_triangle_to_list(f, c, i, abc[2], "c");
	}

	// update struct, and finish
	m->tfirst  = f;
	m->tnext_a = a;
	m->tnext_b = b;
	m->tnext_c = c;
}
int trimesh_get_triangle_fan(int *out, struct trimesh *m, int v)
{
	if (!m->tfirst) trimesh_fill_triangle_fans(m);

	int r = 0;
	int t = m->tfirst[v];
	while (t >= 0) {
		out[r++] = t;
		int T[3];
		sort_three_ints(T, m->t + 3*t);
		if (v == T[0]) t = m->tnext_a[t];
		if (v == T[1]) t = m->tnext_b[t];
		if (v == T[2]) t = m->tnext_c[t];
	}
	return r;
}

static void trimesh_test_triangle_fans(struct trimesh *m)
{
	printf("nv = %d\n", m->nv);
	printf("nt = %d\n", m->nt);
	for (int i = 0; i < m->nt; i++)
		printf("t[%d] = %d %d %d\n",
				i, m->t[3*i+0], m->t[3*i+1], m->t[3*i+2]);
	for (int i = 0; i < m->nv; i++)
		printf("v[%d] = %d\n", i, m->tfirst[i]);
	for (int i = 0; i < m->nt; i++)
		printf("abc[%d] = %d %d %d\n", i,
				m->tnext_a[i], m->tnext_b[i], m->tnext_c[i]);

	for (int i = 0; i < m->nv; i++)
	{
		int out[m->nt];
		int nout = trimesh_get_triangle_fan(out, m, i);
		fprintf(stderr, "trifan[%d] (%d) :", i, nout);
		for (int j = 0; j < nout; j++)
			fprintf(stderr, " %d", out[j]);
		fprintf(stderr, "\n");
	}

}
#endif//TRIMESH_MORE_STUFF


#ifdef TRIMESH_CALCULUS
void trimesh_gradient(                                                   // {{{1
		struct trimesh *m,  // base mesh
		float *y,           // output gradient (array of length m->ne/2)
		float *x            // input function (array of length m->nv)
		)
{
	for (int i = 0; i < m->ne/2; i++)
		y[i] = x[ m->e[2*i+1] ]  -  x[ m->e[2*i+1] ];
}

void trimesh_centering(                                                  // {{{1
		struct trimesh *m,  // base mesh
		float *y,           // output field (array of length m->ne/2)
		float *x            // input function (array of length m->nv)
		)
{
	for (int i = 0; i < m->ne/2; i++)
		y[i] = 0.5 * (  x[ m->e[2*i+1] ]  +  x[ m->e[2*i+1] ]  );
}

void trimesh_recentering(                                                // {{{1
		struct trimesh *m,  // base mesh
		float *y,           // output function (array of length m->nv)
		float *x            // input field (array of length m->ne/2)
		)
{
	for (int i = 0; i < m->ne/2; i++)
		y[i] = 0.5 * (  x[ m->e[2*i+1] ]  +  x[ m->e[2*i+1] ]  );
}

void trimesh_divergence(                                                 // {{{1
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

void trimesh_laplacian(                                                  // {{{1
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


// example "main" function                                                  {{{1
#ifdef TRIMESH_DEMO_MAIN
#include "iio.h"
int main(int c, char *v[])
{
	if (c != 2 && c != 3)
		return fprintf(stderr, "usage:\n\t%s in.tif [out.ply]\n", *v);
	//                                         0 1       2
	char *filename_in  = v[1];
	char *filename_out = c > 2 ? v[2] : "-";

	//// read input DSM
	//int w, h;
	//float *x = iio_read_image_float(filename_in, &w, &h);

	//// create triangulation
	//struct trimesh m[1];
	//trimesh_create_from_dem(m, x, w, h);
	struct trimesh m[1];
	trimesh_read_from_ply(m, filename_in);

#ifdef TRIMESH_MORE_STUFF
	// stuff
	trimesh_fill_triangle_fans(m);
	trimesh_test_triangle_fans(m);
#endif//TRIMESH_MORE_STUFF

	// save triangulation into ply flie
	trimesh_write_to_ply(filename_out, m);

	// cleanup and exit
	trimesh_free_tables(m);
	//free(x);

	// do a silly consistency loop
	if (*filename_out != '-') {
		struct trimesh n[1];
		trimesh_read_from_ply(n, filename_out);
		trimesh_write_to_ply("/tmp/plytmp.ply", n);
		trimesh_free_tables(n);
	}

	// cleanup and exit
	//free(x);
	return 0;
}
#endif//TRIMESH_DEMO_MAIN


// vim:set foldmethod=marker:
