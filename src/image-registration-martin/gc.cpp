#include <iostream>
#include "gradcorr.h"


int main(int argc, char **argv) {
    if (argc != 3) {
        std::cout << "Usage: gc I1 I2" << std::endl;
        std::cout << "I1 and I2 must be grayscale and have the same size." << std::endl << std::endl;
    } else {
        double resX, resY;
        registerGCFiles(argv[1], argv[2], &resX, &resY);
        std::cout << "Estimated displacement: (" << resX << "," << resY << ")" << std::endl << std::endl;
    }
    return 0;
}

