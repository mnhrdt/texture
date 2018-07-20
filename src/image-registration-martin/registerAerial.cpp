#include <iostream>
#include "satreg.h"


int main(int argc, char **argv) {
    if (argc != 4)
        std::cout << "Usage: registerAerial {satImage}.tif {aerialImage}.tif {aerialImage}_PT_opt2ref.txt" << std::endl << std::endl;
    else {
        double resX, resY;
        register_satellite_aerial(argv[1],
                                  argv[2],
                                  argv[3], resX, resY);
        //std::cout << "Res: (" << resX << "," << resY << ")" << std::endl;

    }
    return 0;
}

