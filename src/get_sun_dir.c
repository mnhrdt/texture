#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// compute sun direction from azimuth and elevation
void sun_direction(double n[3], double azimuth, double elevation)
{
    double az_rad = azimuth * M_PI / 180;
    double el_rad = elevation * M_PI / 180;
    n[0] = -cos(el_rad)*sin(az_rad);
    n[1] = -cos(el_rad)*cos(az_rad);
    n[2] = -sin(el_rad);
}

int main(int argc, char *v[])
{
    if (argc < 3)
        return fprintf(stderr, "usage:\n\t"
                "%s azimuth elevation \n",*v);
                //0 1        2
    double  azimuth = atoi(v[1]);
    double  elevation = atoi(v[2]);

    double n[3] = {0, 0, 0};
    sun_direction(n, azimuth, elevation);
    printf("%lf %lf %lf\n", n[0], n[1], n[2]);
    return 0;
}
