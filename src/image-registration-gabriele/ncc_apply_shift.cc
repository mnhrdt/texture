#include <stdio.h>
extern "C" {
    #include "iio.h"
}

// a structure to wrap images
#include "img.h"
#include "img_tools.h"


int main(int argc, char* argv[])
{
    if (argc < 6) {
        fprintf(stderr, "too few parameters\n");
        fprintf(stderr, "\tusage: %s input_img dx dy a b [output_img (stdout)]\n", argv[0]);
        return 1;
    }

    // read the parameters
    char *input_image = argv[1];
    int dx = atoi(argv[2]);
    int dy = atoi(argv[3]);
    float a = atof(argv[4]);
    float b = atof(argv[5]);
    char *output_image = (argc > 6) ? argv[6] : (char *) "-";

    // read input
    struct Img u = iio_read_vector_split(input_image);
    struct Img v(u);

    // apply transformation
    for (int y = 0; y < v.ny; y++)
    for (int x = 0; x < v.nx; x++)
        v(x, y) = a * valnan(u, x+dx, y+dy) + b;

    // write output
    iio_write_vector_split(output_image, v);
    return 0;
}
