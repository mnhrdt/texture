Image registration by gradient phase correlation
================================================

Version 1  (october 11, 2017)
Copyright: Martin Rais 2017



REQUIREMENTS
------------

To compile this code, the following libraries/packages are required: CMake,
libTIFF and FFTW.


COMPILATION
-----------

To generate the binaries, from the folder of CMakeLists.txt, execute the
following:

	mkdir build
	cd build
	cmake ..
	make

This will generate three executables: gc, registerAerial and RegisterAll



USAGE
-----

gc: takes 2 grayscale images of the same size and prints their displacement
(estimated using gradient correlation) on the screen.

registerAerial: takes the satellite image (ortho image), the aerial image and
the text file with the homography and prints the displacement on the screen.

registerAll: takes a single file (input.txt) as parameter and performs all
registrations printing every displacement on the screen and generating an image
with all registered inputs (result.tif).

The file format of input.txt must be the following:
Line 1: absolute path where all images are stored
Line 2: filename (including file extension) of the satellite image (ortho) without the path.
Line 3-N: filenames (including file extensions) of aerial images without the path.


EXAMPLES
--------


$ build/gc data/lena.png data/lena2.png

Image 1 size: 256x256x1
Image 2 size: 256x256x1
Estimated displacement: (1.81511,-3.92628)



$ build/registerAerial data/CMLA_ortho_23_mars_2017.tif data/CMLA_image_075_b2.tif data/CMLA_image_075_b2_PT_opt2ref.txt

Satellite image  size: 1510x1140x3
Aerial image size: 1024x1024
H inv: 4.36728 -1.57096 -1601.98 3.60589 9.03815 -7093.46 0.00014111 -5.60396e-05 1
Estimated displacement: (32.3863,-65.0634)



prints all displacements and generates result.tif with all registered images

$ build/registerAll data/input.txt

