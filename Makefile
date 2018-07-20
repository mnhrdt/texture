F = -march=native -O3 -DNDEBUG
F = -g -Wall -Wno-unused -Werror -Wno-deprecated# -fsanitize=address
A =# -fsanitize=address

# configuration
CFLAGS = $F `gdal-config --cflags` $A
LDLIBS = -lm -ltiff -ljpeg -lpng -lGeographic \
	 -lCGAL -lgmp `gdal-config --libs` $A


# framework compatibility shit (seriously, frameworks are pure brain damage)
ifeq (Darwin,$(shell uname))
  LDLIBS := $(LDLIBS) -framework GLUT -framework OpenGL
else
  LDLIBS := $(LDLIBS) -lglut -lGL
endif


# variables
OBJ = iio drawtriangle trimesh rpc pickopt normals geographiclib_wrapper
BIN = get_utm_normal_shadow colorize_vertices_from_one_image \
      write_coloured_ply refine mflip triproc \
      ncc_compute_shift gc

OBJ := $(OBJ:%=src/%.o)
BIN := $(BIN:%=bin/%)

# default target
all: $(BIN)  # octave_iio

# rule for C sources
bin/% : src/%.c $(OBJ)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

# rule for C++ sources
bin/% : src/%.cc $(OBJ)
	$(CXX) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS) -lboost_system #fml

# compile gabriele's registration
bin/ncc_compute_shift:
	$(MAKE) -C src/image-registration-gabriele -j
	cp src/image-registration-gabriele/ncc_compute_shift $@

# compile martin's registration
bin/gc:
	$(MAKE) -C src/image-registration-martin -j
	cp src/image-registration-martin/gc $@

# octave/matlab iio interface
octave_iio:
	$(MAKE) -C src/octave-iio
	cp src/octave-iio/iio_read.mex* scripts
	cp src/octave-iio/iio_write.mex* scripts

# bureaucracy
clean:
	$(RM) $(OBJ) $(BIN)
	$(RM) scripts/iio_*.mex*
	$(MAKE) -C src/image-registration-martin clean
	$(MAKE) -C src/image-registration-gabriele clean
.SECONDARY : $(OBJ)

# test
test: bin/refine
	./bin/refine data/a.off data/a_out.off --res 0.05

#gltest: ./bin/glflip data/fine_mesh_02.off data/rgb_02.tiff data/rgb_08.tiff
#	$^
