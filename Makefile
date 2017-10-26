# the following option is used to control all C and C++ compilations
FLAGS ?= -g -Wall -Wextra -Werror -Wno-unused

# required libraries
IIOLIBS := -lz -ltiff -lpng -ljpeg -lm
GEOLIBS := -lstdc++ -lGeographic

# variables for implicit rules
CFLAGS = $(FLAGS)
CXXFLAGS = $(FLAGS)
override CPPFLAGS := $(shell gdal-config --cflags)
override LDLIBS   := $(LDLIBS) $(shell gdal-config --libs) $(IIOLIBS) $(GEOLIBS)

# binaries
BIN := colorize get_corners get_P_of_crop colorize_with_shadows colorsingle \
       colormultiple triangles zbuffer get_projection_matrix vector_proj \
       colorfancy get_msi_offset recale
BIN := $(addprefix bin/,$(BIN))
OBJ := src/iio.o src/geographiclib_wrapper.o

# default target (build all the programs inside bin/)
default: $(BIN)

# rule to build all the targets with the same pattern
$(BIN) : bin/% : src/%.o $(OBJ)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS)



# bureaucracy
-include .deps.mk
.deps.mk : ; $(CC) $(CPPFLAGS) -MM src/*.c|sed 's_^[a-z]_src/&_'>$@
clean    : ; $(RM) $(BIN) src/*.o
.PHONY   : clean default


# unit test
test: bin/colorize
	./bin/colorize data/Challenge1_Lidar.tif -21 data/img_01.ntf data/img_01.rpc out_colors.tif


# compatibility hack for older compilers
ifeq (,$(shell $(CC) -dM -E - < /dev/null | grep __STDC_VERSION__))
CFLAGS := $(CFLAGS) -std=gnu99
endif


# other builds
bin/ncc_apply_shift: src/ncc_apply_shift.cc src/img.cc src/point.cc
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $^ src/iio.o -lpng -ltiff -ljpeg -o $@

bin/ncc_compute_shift: src/ncc_compute_shift.cc src/img.cc src/point.cc
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $^ src/iio.o -lpng -ltiff -ljpeg -o $@

bin/trimesh: src/trimesh.c src/iio.o
	$(CC) $(CFLAGS) $(CPPFLAGS) -DTRIMESH_DEMO_MAIN $^ -o $@ $(IIOLIBS)
