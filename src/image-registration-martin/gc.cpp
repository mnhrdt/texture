#include <iostream>
#include <cstring>
#include "gradcorr.h"

/* c pointer to original argc
   v pointer to original argv
   o option name (after hyphen)
   d default value */
static const char *pick_option(int *c, char ***v, const char *o, const char *d)
{
  int argc = *c;
  char **argv = *v;
  int id = d ? 1 : 0;
  int i;
  for (i = 0; i < argc - id; i++)
    if (argv[i][0] == '-' && 0 == strcmp(argv[i]+1, o))
      {
        char *r = argv[i+id]+1-id;
        int j;
        *c -= id+1;
        for (j = i; j < argc - id; j++)
          (*v)[j] = (*v)[j+id+1];
        return r;
      }
  return d;
}

int main(int argc, char **argv) {
    bool unix_style = pick_option(&argc, &argv, "u", "1");
    if (argc != 3) {
        std::cout << "Usage: gc I1 I2" << std::endl;
        std::cout << "I1 and I2 must be grayscale and have the same size." << std::endl << std::endl;
    } else {
        double resX, resY, resZ;
        registerGCFiles(argv[1], argv[2], &resX, &resY);
        height_shift(argv[1], argv[2], resX, resY, &resZ);
	if (unix_style)
		std::cout << resX << " " << resY << " " << -resZ << "\n";
	else
		std::cout << "Estimated displacement: (" << resX << "," << resY << ")" << std::endl << std::endl;
    }
    return 0;
}

