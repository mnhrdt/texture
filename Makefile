#CC=gcc-7 
#CXX=g++-7
# flags
LDLIBS = -lm -ltiff -ljpeg -lpng -lstdc++ -lGeographic `gdal-config --libs`
CFLAGS = `gdal-config --cflags`

# variables
OBJ = iio drawtriangle trimesh rpc pickopt normals geographiclib_wrapper
BIN = get_utm_normal_shadow colorize_vertices_from_one_image write_coloured_ply refine
OBJ := $(OBJ:%=src/%.o)
BIN := $(BIN:%=bin/%)

# default target
all: $(BIN)

# one rule to compile them all
bin/% : src/%.c $(OBJ)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

# bureaucracy
clean: ; $(RM) $(OBJ) $(BIN)

src/refine.o: src/refine.cpp
	$(CXX) -c src/refine.cpp -o src/refine.o -g -Wall -Wextra -Werror -Wno-unused -Wno-unused-parameter #-I/Users/dautume/Documents/doctorat/useful_codes/others/CGAL-4.12/include/ 

bin/refine: src/refine.o
	$(CXX) src/refine.o -lCGAL -o bin/refine -lgmp #-L/Users/dautume/Documents/doctorat/useful_codes/others/CGAL-4.12/lib/

test: bin/refine 
	./bin/refine data/mesh_curve_scaled_remeshed_02.off data/mesh_test_09.off --res 0.9
