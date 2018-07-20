//
// Created by rais on 22/09/17.
//

#ifndef GC_GAUSSFILTER_H
#define GC_GAUSSFILTER_H

#include <cmath>
#include <cstdio>

void filter(double gk[][5])
{
    double stdv = 1.0;
    double r, s = 2.0 * stdv * stdv;  // Assigning standard deviation to 1.0
    double sum = 0.0;   // Initialization of sun for normalization
    for (int x = -2; x <= 2; x++) // Loop to generate 5x5 kernel
    {
        for(int y = -2; y <= 2; y++)
        {
            r = sqrt(x*x + y*y);
            gk[x + 2][y + 2] = (exp(-(r*r)/s))/(M_PI * s);
            sum += gk[x + 2][y + 2];
        }
    }

    for(int i = 0; i < 5; ++i) // Loop to normalize the kernel
        for(int j = 0; j < 5; ++j)
            gk[i][j] /= sum;

}
#endif //GC_GAUSSFILTER_H
