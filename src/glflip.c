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

	// data : a list of 3d points, several lists of RGB colorings
	int npoints;    // number of points in the mesh
	float *point;   // 3 * npoints
	int ncolors;    // number of coloring arrays
	uint8_t *color; // 3 * npoints * ncolors

	// visualization state
	int idx_coloring;   // a number between 0 and ncolorings-1
	float point_radius;

	// trackball
	float view_quat[4];        // quaternion position
	float view_quat_diff[4];   // quaternion speed
	float view_scale;
	float view_center[3];
	float view_distance;
	float trackball_begin_x;
	float trackball_begin_y;
	float trackball_dx;
	float trackball_dy;
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

static void draw_pped(GLfloat *a, GLfloat *b, bool solid)
{
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
		glNormal3fv(n[i]);
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
	glClear(GL_COLOR_BUFFER_BIT);

	glLoadIdentity();
	glTranslatef(0, 0, e->view_distance);
	glScalef(e->view_scale, e->view_scale, e->view_scale);

	float m[4][4];
	add_quats (e->view_quat_diff, e->view_quat, e->view_quat);
	build_rotmatrix (m, e->view_quat);
	glMultMatrixf (&m[0][0]);

	glBegin(GL_POINTS);
	for (int i = 0; i < e->npoints; i++)
	{
		uint8_t *c = (e->color + 3*e->npoints*e->idx_coloring) + 3*i;
		float   *p = e->point + 3*i;
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


//	GLfloat a[3] = {0, 0, 0};
//	GLfloat b[3] = {1, 1, 1};
//	draw_pped(a, b, false);

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
		e->view_quat_diff[0] = 0;
		e->view_quat_diff[1] = 0;
		e->view_quat_diff[2] = 0;
		e->view_quat_diff[3] = 1;
		e->trackball_begin_x = x;
		e->trackball_begin_y = y;
	}

	if (b == FTR_BUTTON_LEFT && s==1)
	{
		e->view_quat_diff[0] = 0;
		e->view_quat_diff[1] = 0;
		e->view_quat_diff[2] = 0;
		e->view_quat_diff[3] = 1;
		e->trackball_dx = 0.0;
		e->trackball_dy = 0.0;
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

	trackball(e->view_quat,
			(2.0 * e->trackball_begin_x - w) / w,
			(h - 2.0 * e->trackball_begin_y) / h,
			(2.0 * x - w) / w,
			(h - 2.0 * y) / h);

	e->trackball_dx = x - e->trackball_begin_x;
	e->trackball_dy = y - e->trackball_begin_y;
	glutPostRedisplay();
}

// glut handler : glutPassiveMotionFunc
void my_passivemotionfunc(int x, int y)
{
	fprintf(stderr, "GLUT passive motion %d %d\n", x, y);
}


// unified key handler for special and non-special keys
static void key_handler(int k, int x, int y)
{
	struct state *e = &global_state;

	if (k == FTR_KEY_ESC || k == 'q')
		glutLeaveMainLoop();

	if (k == 'd') e->view_distance += 0.1;
	if (k == 'D') e->view_distance -= 0.1;
	if (k == 's') e->view_scale *= 1.3;
	if (k == 'S') e->view_scale /= 1.3;
	if (k == ' ') e->idx_coloring = (1 + e->idx_coloring) % e->ncolors;

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
	//glutReshapeFunc(my_reshapefunc);
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
}


// silly function to create a cube with colored sides
static void fill_synthetic_rgbcube(struct state *e,
		int n,  // number of points per side
		float s // length of the side
		)
{
	float z = 255;                          // top of color range
	float   *p = malloc(12*3*n*sizeof*p);   // array of 3D points
	uint8_t *c = malloc(2*12*3*n*sizeof*c); // array RGB colors
	float d = s / (n - 1);                  // step in 3D space
	float k = z / (n - 1);                  // step in RGB space

	// todo: break the loop into several separate loops for better ordering
	for (int i = 0; i < n; i++)
	{
		fprintf(stderr, "i = %d, p = %p\n", i, (void*)p);
		// (i*d, -, -) : (i*c, -, -)
		*p++ = i*d; *p++ = 0; *p++ = 0; *c++ = i*k; *c++ = 0; *c++ = 0;
		*p++ = i*d; *p++ = s; *p++ = 0; *c++ = i*k; *c++ = z; *c++ = 0;
		*p++ = i*d; *p++ = 0; *p++ = s; *c++ = i*k; *c++ = 0; *c++ = z;
		*p++ = i*d; *p++ = s; *p++ = s; *c++ = i*k; *c++ = z; *c++ = z;

		// (-, i*d, -) : (-, i*c, -)
		*p++ = 0; *p++ = i*d; *p++ = 0; *c++ = 0; *c++ = i*k; *c++ = 0;
		*p++ = s; *p++ = i*d; *p++ = 0; *c++ = z; *c++ = i*k; *c++ = 0;
		*p++ = 0; *p++ = i*d; *p++ = s; *c++ = 0; *c++ = i*k; *c++ = z;
		*p++ = s; *p++ = i*d; *p++ = s; *c++ = z; *c++ = i*k; *c++ = z;

		// (-, -, i*d) : (-, -, i*c)
		*p++ = 0; *p++ = 0; *p++ = i*d; *c++ = 0; *c++ = 0; *c++ = i*k;
		*p++ = s; *p++ = 0; *p++ = i*d; *c++ = z; *c++ = 0; *c++ = i*k;
		*p++ = 0; *p++ = s; *p++ = i*d; *c++ = 0; *c++ = z; *c++ = i*k;
		*p++ = s; *p++ = s; *p++ = i*d; *c++ = z; *c++ = z; *c++ = i*k;
	}
	for (int i = 0; i < 12*n; i++)
	{
		*c++ = 255;
		*c++ = 0;
		*c++ = 0;
	}

	e->point = p - 12*3*n;
	e->color = c - 2*12*3*n;
	e->npoints = 12*n;
	e->ncolors = 2;
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
	e->view_scale = 1;
	e->view_distance = 0;

	e->point_radius = 1;

	fill_synthetic_rgbcube(e, 60, 1);
}

int main(void)
{
	struct state *e = &global_state;
	setup_initial_state(e, 800, 600);
	setup_glut_environment(e->w, e->h);
	glutMainLoop();
	return 0;
}
