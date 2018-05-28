# flags
F = -g -fsanitize=address -Wall -Wno-unused -Werror
F = -march=native -O3 -DNDEBUG 

# configuration
CFLAGS = $F `gdal-config --cflags`
LDLIBS = -lm -ltiff -ljpeg -lpng -lGeographic \
	 -lCGAL -lgmp `gdal-config --libs` $F

# variables
OBJ = iio drawtriangle trimesh rpc pickopt normals geographiclib_wrapper
BIN = get_utm_normal_shadow colorize_vertices_from_one_image \
      write_coloured_ply refine

OBJ := $(OBJ:%=src/%.o)
BIN := $(BIN:%=bin/%)

# default target
all: $(BIN)

# rule for C sources
bin/% : src/%.c $(OBJ)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

# rule for C++ sources
bin/% : src/%.cpp $(OBJ)
	$(CXX) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS) -lboost_system #fml

# bureaucracy
clean: ; $(RM) $(OBJ) $(BIN)

# test
test: bin/refine
	./bin/refine data/a.off data/a_out.off --res 0.05
