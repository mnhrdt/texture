#include <iostream>
#include <memory.h>
#include <vector>
#include <cmath>
#include "satreg.h"
extern "C" {
#include "iio.h"
}

int main(int argc, char **argv) {
    if (argc != 2) {
        std::cout << "Usage: convertAll input.txt" << std::endl << std::endl;
        std::cout <<
                  "input.txt file format should be:"
                  << std::endl
                  << "Line 1: directory where all filenames are stored" << std::endl
                  << "Line 2: filename of satellite ortho image" << std::endl
                  << "Line 3-N: filenames of aereal images" << std::endl;
    }
    else {
        FILE *f = fopen(argv[1], "rt");
        char line[1024];
        char *token;
        int i = 0;
        char satImgFile[100];
        char aerialFile[100];
        char path[1024];
        char tmpFilename[1024];
        double *finalImage = NULL;
        std::vector<double *> registeredImages;
        int w, h;
        double resX, resY;
        while (fgets(line, 1024, f))
        {
            if (line[0] !=  '\n') {
                if (i==0) {
                    strcpy(path, line);
                    strtok(path, "\n");
                } else  if (i==1) {
                    strcpy(satImgFile, line);
                    strtok(satImgFile, "\n");
                    strcpy(tmpFilename, path);
                    strcat(tmpFilename, satImgFile);
                    std::cout << "Processing " << tmpFilename << std::endl;
                    // Open sat image to initialize result
                    int pd;
                    double *img = iio_read_image_double_vec(tmpFilename, &w, &h, &pd);
                    free(img);
                } else {
                    strcpy(aerialFile, line);
                    strtok(aerialFile, "\n");
                    strcpy(tmpFilename, path);
                    strcat(tmpFilename, aerialFile);

                    std::cout << "Processing " << tmpFilename << std::endl;
                    char aerialGTFile[1024];
                    const char *last = strrchr(tmpFilename, '.');
                    int length = last - tmpFilename;
                    strncpy(aerialGTFile, tmpFilename, length);
                    aerialGTFile[length] = '\x0';
                    strcat(aerialGTFile, "_PT_opt2ref.txt");
                    std::cout << aerialGTFile << std::endl;
                    register_satellite_aerial(satImgFile,
                                              aerialFile,
                                              aerialGTFile, resX, resY,
                                              &finalImage);
                    registeredImages.push_back(finalImage);
                }
            }
            i++;
        }
        fclose(f);
        std::cout << "Merging results into single image..." << std::endl;
        // Merge results into a single image
        double *fusedImage = (double *) calloc(w * h, sizeof(double));
        double *acumMat = (double *) calloc(w * h, sizeof(double));
        for (int i=0;i<registeredImages.size();i++) {
            double *img = registeredImages[i];
            for (int x = 0; x < w; x++) {
                for (int y = 0; y < h; y++) {
                    int pos = x + w * y;
                    double val = img[pos];
                    if (!std::isnan(val)) {
                        fusedImage[pos] += val;
                        acumMat[pos]++;
                    }
                }
            }
        }
        for (int x=0;x<w;x++) {
            for (int y = 0; y < h; y++) {
                if (acumMat[x+y*w] != 0)
                    fusedImage[x + y*w] /= acumMat[x+y*w];
            }
        }
        iio_write_image_double("result.tif", fusedImage, w, h);
        std::cout << "Registration finished correctly." << std::endl;
        for (auto img : registeredImages)
            free(img);


    }
    return 0;
}
