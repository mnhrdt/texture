// silly example for using the "traverse_triangle" function
//
// fill an image with triangles at various positions
// there should be no holes betwen the triangles

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include "drawtriangle.c"

struct image {
	int w, h;    // image width, height
	uint8_t *x;  // image data
	int r, g, b; // current color
};

// put a pixel at the specified position
// (ee must point to an allocated image)
static void putpixel(int i, int j, void *ee)
{
	struct image *e = ee;
	if (i >= 0 && i < e->w && j >= 0 && j < e->h)
	{
		int idx = j*e->w + i;
		e->x[3*idx+0] = e->r;
		e->x[3*idx+1] = e->g;
		e->x[3*idx+2] = e->b;
	}
}

// return a random number between 0 and 1
static float random_uniform(void)
{
	return rand()/(0.0+RAND_MAX);
}

#include "iio.h"
int main(void)
{
	// create a black image
	int w = 512, h = 512;
	uint8_t *x = malloc(3*w * h * sizeof*x);
	for (int i = 0; i < 3*w*h; i++)
		x[i] = 0;

	// create a grid of points
	int nx = 10;
	int ny = 10;
	float grid[ny*nx][2];
	for (int j = 0; j < ny; j++)
	for (int i = 0; i < nx; i++)
	{
		grid[j*nx+i][0] = 30*random_uniform() + i * w / (float)nx;
		grid[j*nx+i][1] = 30*random_uniform() + j * h / (float)ny;
	}

	// join these points by triangles
	int nt = 2 * (nx-1) * (ny-1), cx = 0;
	int tri[nt][3];
	for (int j = 0; j < ny-1; j++)
	for (int i = 0; i < nx-1; i++)
	{
		tri[cx][0] = (j+0)*nx+i+0;
		tri[cx][1] = (j+0)*nx+i+1;
		tri[cx][2] = (j+1)*nx+i+0;
		cx += 1;
		tri[cx][0] = (j+1)*nx+i+1;
		tri[cx][1] = (j+0)*nx+i+1;
		tri[cx][2] = (j+1)*nx+i+0;
		cx += 1;
	}
	assert(cx == nt);

	// fill-in the image
	for (int i = 0; i < nt; i++)
	{
		float *a = grid[tri[i][0]];
		float *b = grid[tri[i][1]];
		float *c = grid[tri[i][2]];
		float abc[3][2] = { {a[0], a[1]}, {b[0], b[1]}, {c[0], c[1]} };
		struct image e = {
			.w = w,
			.h = h,
			.x = x,
			.r = 255*random_uniform(),
			.g = 255*random_uniform(),
			.b = 255*random_uniform(),
		};
		traverse_triangle(abc, putpixel, &e);
	}
	iio_write_image_uint8_vec("triangles.png", x, w, h, 3);
	fprintf(stderr, "created images \"triangles.png\"\n");
	free(x);
	return 0;
}
