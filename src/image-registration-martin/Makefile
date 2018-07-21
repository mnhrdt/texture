# choice between debug and relsease (the second line prevails)
FLAGS = -g -Wall -Wextra -fsanitize=address            # debug
#FLAGS = -O3 -march=native -DNDEBUG  # release

# compiler flags
CXXFLAGS ?= $(FLAGS)
CFLAGS   ?= $(FLAGS)

# required libraries
LDLIBS = -ltiff -ljpeg -lpng -lfftw3 -lz -lm

# variables
BIN    = registerAerial gc
OBJ    = gradcorr.o gradient.o iio.o interpolation.o matops.o satreg.o window.o

# default target (all the binaries)
all    : $(BIN)

# dependences (each binary depends on all the objects)
$(BIN) : $(OBJ)

# bureaucracy
clean  :
	$(RM) $(BIN) $(OBJ)
