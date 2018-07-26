#include "trimesh.h"

#include <stdio.h>
#include <string.h>
int main_off2edges(int c, char *v[])
{
	if (c != 3)
		return fprintf(stderr, "usage:\t%s in.off out.txt\n", *v);
	//                                       0 1      2
	char *filename_in  = v[1];
	char *filename_out = v[2];

	struct trimesh m[1];
	trimesh_read_from_off(m, filename_in);
	trimesh_dump_edges(filename_out, m);
	trimesh_free_tables(m);

	return 0;
}

int main(int c, char *v[])
{
	if (c < 2)
		return fprintf(stderr, "usage:\t%s {off2edges|...}\n", *v);
	if (0 == strcmp(v[1], "off2edges")) return main_off2edges(c-1, v+1);
	return 1;
}
