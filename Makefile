# the following two options are used to control all C and C++ compilations
CFLAGS   ?= -march=native -O3
CXXFLAGS ?= -march=native -O3

# required libraries
IIOLIBS := -lz -ltiff -lpng -ljpeg -lm
GEOLIBS := -lstdc++ -lGeographic

# variables
override CFLAGS := $(CFLAGS) `gdal-config --cflags`
override LDLIBS := $(LDLIBS) `gdal-config --libs` $(IIOLIBS) $(GEOLIBS)

# binaries
BIN := colorize get_corners get_P_of_crop colorize_with_shadows colorsingle triangles
BIN := $(addprefix bin/,$(BIN))
OBJ := src/iio.o src/geographiclib_wrapper.o

# default target (build all the programs inside bin/)
default: $(BIN)

# rule to build all the targets with the same pattern
$(BIN) : bin/% : src/%.o $(OBJ)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS)

# bureaucracy
clean: ; $(RM) $(BIN) src/*.o
.PHONY: clean default

# unit test
test: bin/colorize
	./bin/colorize data/Challenge1_Lidar.tif -21 data/img_01.ntf data/img_01.rpc out_colors.tif

# compatibility hack for older compilers
ifeq (,$(shell $(CC) -dM -E - < /dev/null | grep __STDC_VERSION__))
CFLAGS := $(CFLAGS) -std=gnu99
endif

