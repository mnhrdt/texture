#ifndef PICKOPT_H
#define PICKOPT_H

#include <string.h>
#include <stdbool.h>

// @c pointer to original argc
// @v pointer to original argv
// @o option name (after hyphen)
// @d default value
char * pick_option(int *c, char ***v, char *o, char *d);

#endif//PICKOPT_H
