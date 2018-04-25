#ifndef DRAWTRIANGLE_H
#define DRAWTRIANGLE_H

// fill a triangle defined by three points a,b,c
void traverse_triangle(
		float abc[3][2],          // coordinates of the three points
		void (*f)(int,int,void*), // function to apply to each pixel
		void *e                   // passed to f
		);

#endif//DRAWTRIANGLE_H
