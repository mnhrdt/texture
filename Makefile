CC=gcc-7 
CXX=g++-7
# flags
LDLIBS = -lm -ltiff -ljpeg -lpng -lstdc++ -lGeographic `gdal-config --libs`
CFLAGS = `gdal-config --cflags`

# variables
OBJ = iio drawtriangle trimesh rpc pickopt normals geographiclib_wrapper
BIN = get_utm_normal_shadow colorize_vertices_from_one_image write_coloured_ply fusion
OBJ := $(OBJ:%=src/%.o)
BIN := $(BIN:%=bin/%)

# default target
all: $(BIN)

# one rule to compile them all
bin/% : src/%.c $(OBJ)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

# bureaucracy
clean: ; $(RM) $(OBJ) $(BIN)

test: bin/refine data/mesh_curve_scaled_remeshed_02.off data/mesh_test_09.off --res 0.9
