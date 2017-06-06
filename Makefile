# the following two options are used to control all C and C++ compilations
export CFLAGS =   -march=native -O3
export CXXFLAGS = -march=native -O3

# these options are only used for the programs directly inside "./c/"
IIOLIBS = -lz -ltiff -lpng -ljpeg -lm

# The following conditional statement appends "-std=gnu99" to CFLAGS when the
# compiler does not define __STDC_VERSION__.  The idea is that many older
# compilers are able to compile standard C when given that option.
# This hack seems to work for all versions of gcc, clang and icc.
CVERSION = $(shell $(CC) -dM -E - < /dev/null | grep __STDC_VERSION__)
ifeq ($(CVERSION),)
CFLAGS := $(CFLAGS) -std=gnu99
endif

# default rule
default: bin/colorize bin/get_corners bin/get_P_of_crop bin/colorize_with_shadows

# rules to build each program

bin/colorize: src/iio.o src/geographiclib_wrapper.o src/colorize.c
	$(CC) $(CFLAGS) `gdal-config --cflags` $^ $(IIOLIBS) -lstdc++ -lGeographic `gdal-config --libs` -o $@

bin/colorize_with_shadows: src/iio.o src/geographiclib_wrapper.o src/colorize_with_shadows.c
	$(CC) $(CFLAGS) `gdal-config --cflags` $^ $(IIOLIBS) -lstdc++ -lGeographic `gdal-config --libs` -o $@

bin/get_corners: src/iio.o src/geographiclib_wrapper.o src/get_corners.c
	$(CC) $(CFLAGS) `gdal-config --cflags` $^ $(IIOLIBS) -lstdc++ -lGeographic `gdal-config --libs` -o $@

bin/get_P_of_crop: src/iio.o src/geographiclib_wrapper.o src/get_P_of_crop.c
	$(CC) $(CFLAGS) `gdal-config --cflags` $^ $(IIOLIBS) -lstdc++ -lGeographic `gdal-config --libs` -o $@

clean:
	$(RM) src/*.o bin/colorize out_colors.tif

test: bin/colorize
	./bin/colorize data/Challenge1_Lidar.tif -21 data/img_01.ntf data/img_01.rpc out_colors.tif
