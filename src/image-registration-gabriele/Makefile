CXXFLAGS = -fopenmp -O3 -Wall -Wextra -Wno-unused -Wno-unused-parameter
LDLIBS = -lz -ltiff -lpng -ljpeg -lm

BIN = ncc_apply_shift ncc_compute_shift
OBJ = iio.o img.o point.o

default: $(BIN)

$(BIN) : $(OBJ)



## bureaucracy
#-include .deps.mk
#.deps.mk : ; $(CC) $(CPPFLAGS) -MM src/*.c|sed 's_^[a-z]_src/&_'>$@
clean    : ; $(RM) $(BIN) $(OBJ)
.PHONY   : clean default
