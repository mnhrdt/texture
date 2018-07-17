#include <ctype.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <GL/freeglut.h>

// state of the viewer
struct state {
	// window
	int w; // window width;
	int h; // window height;

	// data : a list of 3d points, several lists of RGB colorings, triangles
	int npoints;    // number of points in the mesh
	int ncolors;    // number of coloring arrays
	int ntriangles; // number of triangles in the mesh
	float *point;   // 3 * npoints
	float *color;   // 3 * npoints * ncolors
	int *triangle;  // 3 * ntriangles
	float bbx[6];   // (rectangular bounding box, computed from data)

	// visualization state
	int toggle_triangles; // whether to show triangles or not
	int idx_coloring;   // a number between 0 and ncolorings-1
	float point_radius;
	float center[3];

	// trackball
	float view_quat[4];        // quaternion position
	float view_quat_diff[4];   // quaternion speed
	float view_scale;
	float view_center[3];
	float view_distance;
	float trackball_begin_x;
	float trackball_begin_y;
	float track_motion_x;
	float track_motion_y;
	float track_prev_x;
	float track_prev_y;
	int track_drag;
	float view_cnt, view_cnt_diff; // debugging
};

static struct state global_state; // hack


// regular keys
#define FTR_KEY_ESC        27
#define FTR_KEY_DEL        127

// special keys
#define FTR_KEY_LEFT       100
#define FTR_KEY_UP         101
#define FTR_KEY_RIGHT      102
#define FTR_KEY_DOWN       103
#define FTR_KEY_PAGE_UP    104
#define FTR_KEY_PAGE_DOWN  105
#define FTR_KEY_HOME       106
#define FTR_KEY_END        107
#define FTR_KEY_INSERT     108

// button identifiers
#define FTR_BUTTON_LEFT    256
#define FTR_BUTTON_MIDDLE  512
#define FTR_BUTTON_RIGHT   1024
#define FTR_BUTTON_UP      2048
#define FTR_BUTTON_DOWN    4096


#include "trackball.c" // quaternion rotations and related utilities

static void draw_pped(float *a, float *b, bool solid)
{
	fprintf(stderr, "pped = %g %g %g  %g %g %g\n",
			a[0], a[1], a[2], b[0], b[1], b[2]);
	GLenum type = solid ? GL_QUADS : GL_LINE_LOOP;
	GLfloat v[8][3] =
	{
		{a[0], a[1], a[2]}, // 0
		{b[0], a[1], a[2]}, // 1
		{a[0], b[1], a[2]}, // 2
		{b[0], b[1], a[2]}, // 3
		{a[0], a[1], b[2]}, // 4
		{b[0], a[1], b[2]}, // 5
		{a[0], b[1], b[2]}, // 6
		{b[0], b[1], b[2]}  // 7
	};
	GLint f[6][4] =
	{
		{0, 2, 3, 1},
		{1, 3, 7, 5},
		{4, 5, 7, 6},
		{0, 4, 6, 2},
		{0, 1, 5, 4},
		{2, 6, 7, 3}
	};
	GLfloat n[6][3] =
	{
		{0, 0, -1},
		{1, 0, 0},
		{0, 0, 1},
		{-1, 0, 0},
		{0, -1, 0},
		{0, 1, 0}
	};
	for (int i = 0; i < 6; i++)
	{
		glBegin(type);
		if (solid) glNormal3fv(n[i]);
		for (int j = 0; j < 4; j++)
			glVertex3fv(v[f[i][j]]);
		glEnd();
	}
}


// glut handler : glutDisplayFunc
void my_displayfunc(void)
{
	// TODO: change this to use the modern rendering pipeline


	struct state *e = &global_state;
	fprintf(stderr, "GLUT display\n");

	fprintf(stderr, "\tscale=%g distance=%g\n", e->view_scale, e->view_distance);
	fprintf(stderr, "\tQ = [ %g %g %g %g ; %g %g %g %g ]\n",
			e->view_quat[0], e->view_quat[1],
			e->view_quat[2], e->view_quat[3],
			e->view_quat_diff[0], e->view_quat_diff[1],
			e->view_quat_diff[2], e->view_quat_diff[3]);


	glClearColor(0.5, 0.8, 0.5, 1); // light greenish
	//glClearColor(0, 0, 0, 1); // black
	glClearDepth(1);
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT|GL_STENCIL_BUFFER_BIT);

	glLoadIdentity();

	// display control points in image space
	if (e->track_drag) {
		glPointSize( 7 );
		glBegin(GL_POINTS);
		glColor3ub(255, 0, 0); // red
		glVertex3f(
				(2.0*e->trackball_begin_x - e->w) / e->w,
				(e->h - 2.0*e->trackball_begin_y) / e->h,
				0);
		glColor3ub(255, 255, 255); // white
		glVertex3f(
				(2.0*e->track_prev_x - e->w) / e->w,
				(e->h - 2.0*e->track_prev_y) / e->h,
				0);
		glColor3ub(255, 0, 255); // magenta
		glVertex3f(
				(2.0*e->track_motion_x - e->w) / e->w,
				(e->h - 2.0*e->track_motion_y) / e->h,
				0);
		glEnd();
	}

	// display rest of the points in image space
	glTranslatef(0, 0, e->view_distance);
	glScalef(e->view_scale, e->view_scale, e->view_scale);
	glTranslatef(-e->center[0], -e->center[1], e->center[2]);

	float m[4][4];
	build_rotmatrix (m, e->view_quat);
	glMultMatrixf (&m[0][0]);

	glEnable(GL_ALPHA_TEST);
	glAlphaFunc(GL_NOTEQUAL, 0);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_POINT_SMOOTH);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);

	glPointSize( e->point_radius );
	glBegin(GL_POINTS);
	for (int i = 0; i < e->npoints; i++)
	{
		float *c = (e->color + 3*e->npoints*e->idx_coloring) + 3*i;
		float *p = e->point + 3*i;
		glColor3ub(c[0], c[1], c[2]);
		glVertex3f(p[0], p[1], p[2]);
	}
	glEnd();

//	glBegin(GL_POLYGON);
//	glVertex3f(0.0, 0.0, 0.0);
//	glVertex3f(0.5, 0.0, 0.0);
//	glVertex3f(0.5, 0.5, 0.0);
//	glVertex3f(0.0, 0.5, 0.0);
//	glEnd();


	float cc[6] = {
		e->center[0]-1, e->center[1]-1, e->center[2]-1,
		e->center[0]+1, e->center[1]+1, e->center[2]+1};
	draw_pped(e->bbx, e->bbx+3, false);
	draw_pped(cc, cc + 3, false);

	glFlush();

	glutSwapBuffers();
}


static int button_from_glut_to_ftr(int b)
{
	return 1 << ( b + 8 );
}


// glut handler : glutMouseFunc
void my_mousefunc(int b, int s, int x, int y)
{
	struct state *e = &global_state;
	fprintf(stderr, "GLUT mouse %d %d (%d %d)\n", b, s, x, y);
	b = button_from_glut_to_ftr(b);

	if (b == FTR_BUTTON_LEFT && s==0)
	{
		e->track_drag = 1;
		e->trackball_begin_x = x;
		e->trackball_begin_y = y;
		e->track_motion_x = x;
		e->track_motion_y = y;
	}

	if (b == FTR_BUTTON_LEFT && s==1)
	{
		e->track_drag = 0;
	}
}

// glut handler : glutMotionFunc
void my_motionfunc(int x, int y)
{
	static float ox = 0;
	static float oy = 0;

	fprintf(stderr, "GLUT motion %d %d\n", x, y);

	struct state *e = &global_state;
	float w = e->w;
	float h = e->h;

	e->track_prev_x = e->track_motion_x;
	e->track_prev_y = e->track_motion_y;
	e->track_motion_x = x;
	e->track_motion_y = y;

	fprintf(stderr, "prev=%g %g    mot=%g %g\n",
			e->track_prev_x,
			e->track_prev_y,
			e->track_motion_x,
			e->track_motion_y);

	trackball(e->view_quat_diff,
			(2.0 * e->track_prev_x - w) / w,
			(h - 2.0 * e->track_prev_y) / h,
			(2.0 * x - w) / w,
			(h - 2.0 * y) / h);
        e->trackball_begin_x = x;
        e->trackball_begin_y = y;
	add_quats(e->view_quat, e->view_quat, e->view_quat_diff);

	glutPostRedisplay();
}

// glut handler : glutPassiveMotionFunc
void my_passivemotionfunc(int x, int y)
{
	fprintf(stderr, "GLUT passive motion %d %d\n", x, y);
}

static void action_distance_add(struct state *e, float d)
{
	e->view_distance += d;
}

static void action_scale_mul(struct state *e, float f)
{
	e->view_scale += f;
}

static void action_coloring_cycle(struct state *e)
{
	e->idx_coloring = (1 + e->idx_coloring) % e->ncolors;
}

static void action_point_radius_mul(struct state *e, float f)
{
	e->point_radius *= f;
}


// unified key handler for special and non-special keys
static void key_handler(int k, int x, int y)
{
	struct state *e = &global_state;

	if (k == FTR_KEY_ESC || k == 'q')
		glutLeaveMainLoop();

	if (k == 'd') action_distance_add(e, 0.1);
	if (k == 'D') action_distance_add(e, -0.1);
	if (k == 's') action_scale_mul(e, 1.3);
	if (k == 'S') action_scale_mul(e, 1/1.3);
	if (k == ' ') action_coloring_cycle(e);
	if (k == 'p') action_point_radius_mul(e, 1.1);
	if (k == 'P') action_point_radius_mul(e, 1/1.1);

	glutPostRedisplay();
}

// glut handler : glutKeyboardFunc
void my_keyboardfunc(unsigned char k, int x, int y)
{
	fprintf(stderr, "GLUT keyboard  %d '%c' %d %d\n",
			k, isalnum(k)?k:' ', x, y);

	key_handler(k, x, y);
}

// glut handler : glutSpecialFunc
void my_specialfunc(int k, int x, int y)
{
	fprintf(stderr, "GLUT special  %d '%c' %d %d\n",
			k, isalnum(k)?k:' ', x, y);

	key_handler(1000+k, x, y);
}


void my_reshapefunc(int w, int h)
{
	fprintf(stderr, "GLUT reshape %d %d\n", w, h);

	glViewport (0, 0, w, h);

	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();
	glFrustum(-1, 1,  -1, 1,  1, 60);
	//if (w > h)
	//{
	//	float aspect = w / (float)h;
	//	glFrustum (-aspect, aspect, -1.0, 1.0, 5.0, 60.0);
	//}
	//else
	//{
	//	float aspect = h / (float)w;
	//	glFrustum (-1.0, 1.0, -aspect, aspect, 5.0, 60.0);
	//}

	glMatrixMode (GL_MODELVIEW);
	int e = glGetError();
	if (e)
		fprintf(stderr, "ERR = %d\n", e);
}

// function to call the shitty glut initialization code
static void setup_glut_environment(int w, int h)
{
	int argc = 2;
	char *argv[]={"dummy", "-gldebug", NULL};
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGB|GLUT_DEPTH);
	glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,
			GLUT_ACTION_GLUTMAINLOOP_RETURNS);
	glutCreateWindow("Whatever");
	glutReshapeWindow(w, h);
	glutDisplayFunc(my_displayfunc);
	glutReshapeFunc(my_reshapefunc);
	glutIdleFunc(NULL);
	glutMouseFunc(my_mousefunc);
	glutMotionFunc(my_motionfunc);
	glutPassiveMotionFunc(my_passivemotionfunc);
	glutKeyboardFunc(my_keyboardfunc);
	glutSpecialFunc(my_specialfunc);
	//f->glut_initted = 1;
	//f->glut_button_mask = 0; // not sure
	//f->glut_keymod_mask = 0;
	//fprintf(stderr, "setup glut rgb = %p\n", f->rgb);
	//if (f->glut_window_x >= 0)
	//	glutPositionWindow(f->glut_window_x, f->glut_window_y);
	//
}



static void setup_initial_state(struct state *e, int w, int h)
{
	e->w = w;
	e->h = h;

	e->view_cnt = e->view_cnt_diff = 0;

	e->view_quat[0] = e->view_quat_diff[0] = 0.0;
	e->view_quat[1] = e->view_quat_diff[1] = 0.0;
	e->view_quat[2] = e->view_quat_diff[2] = 0.0;
	e->view_quat[3] = e->view_quat_diff[3] = 1.0;
	e->view_scale = 3;
	e->view_distance = 0;

	e->trackball_begin_x = e->trackball_begin_y = NAN;
	e->track_motion_x = e->track_motion_y = NAN;
	e->track_prev_x = e->track_prev_y = NAN;
	e->track_drag = 0;

	e->point_radius = 1;

	//fill_synthetic_rgbcube(e, 60, 1);
	//fill_synthetic_surface(e, 200, 1);
	e->center[0] = (e->bbx[0] + e->bbx[3])/2;
	e->center[1] = (e->bbx[1] + e->bbx[4])/2;
	e->center[2] = e->bbx[2];
}

#include "iio.h"
#include "trimesh.h"
static void load_data_into_state(struct state *e, char *fname_m,
		int ncolorfiles, char **fnames_c)
{
	// read mesh data
	struct trimesh m[1];
	trimesh_read_from_off(m, fname_m);

	// read color data
	int ncolors = 0;
	int h[ncolorfiles];
	float *tmp[ncolorfiles];
	for (int k = 0; k < ncolorfiles; k++)
	{
		int w, pd;
		tmp[k] = iio_read_image_float_vec(fnames_c[k], &w, h + k, &pd);
		if (pd != 3 || w != m->nv)
			exit(fprintf(stderr, "image \"%s\" dimensions do not "
			       " agree with mesh img=%dx%dx%d, m=%d\n",
			       fnames_c[k], w, h[k], pd, m->nv));
		ncolors += h[k];
	}

	// fill-in data into surface
	e->npoints = m->nv;
	e->ncolors = ncolors;
	e->ntriangles = m->nt;
	e->point = malloc(3*e->npoints * sizeof*e->point);
	e->color = malloc(3*e->npoints * e->ncolors * sizeof*e->color);
	e->triangle = malloc(3*e->ntriangles * sizeof*e->triangle);
	for (int i = 0; i < 3*e->npoints; i++)
		e->point[i] = m->v[i];
	for (int i = 0; i < 3*e->ntriangles; i++)
		e->triangle[i] = m->t[i];
	float *p = e->color;
	for (int k = 0; k < ncolorfiles; k++)
	for (int i = 0; i < 3*e->npoints*h[k]; i++)
		*p++ = tmp[k][i];

	// cleanup
	for (int k = 0; k < ncolorfiles; k++)
		free(tmp[k]);
	trimesh_free_tables(m);

	// compute bounding box
	e->bbx[0] = e->bbx[1] = e->bbx[2] = INFINITY;
	e->bbx[3] = e->bbx[4] = e->bbx[5] = -INFINITY;
	for (int i = 0; i < e->npoints; i++)
	{
		float *x = e->point + 3*i;
		for (int j = 0; j < 3; j++)
		{
			e->bbx[j+0] = fmin(e->bbx[j+0], x[j]);
			e->bbx[j+3] = fmax(e->bbx[j+3], x[j]);
		}
	}
}

int main(int c, char *v[])
{
	// extract input arguments
	if (c < 3)
		return fprintf(stderr, "usage:\n\t%s mesh color1 ...\n", *v);
	//                                         0 1    2 
	char *fname_mesh = v[1];
	int ncolorfiles = c - 2;
	char **fnames_color = v + 2;

	// initialize state
	struct state *e = &global_state;
	load_data_into_state(e, fname_mesh, ncolorfiles, fnames_color);
	setup_initial_state(e, 512, 512);

	// run the viewer
	setup_glut_environment(e->w, e->h);
	glutMainLoop();
	return 0;
}
