F = -march=native -O3 -DNDEBUG
F = -g -Wall -Wno-unused -Werror# -fsanitize=address
A =# -fsanitize=address

# configuration
CFLAGS = $F `gdal-config --cflags` $A
LDLIBS = -lm -ltiff -ljpeg -lpng -lGeographic \
	 -lCGAL -lgmp `gdal-config --libs` -lglut -framework OpenGL $A

# variables
OBJ = iio drawtriangle trimesh rpc pickopt normals geographiclib_wrapper
BIN = get_utm_normal_shadow colorize_vertices_from_one_image \
      write_coloured_ply refine glflip mflip

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
.SECONDARY : $(OBJ)

# test
test: bin/refine
	./bin/refine data/a.off data/a_out.off --res 0.05

gltest: ./bin/glflip data/fine_mesh_02.off data/rgb_02.tiff data/rgb_08.tiff
	$^
