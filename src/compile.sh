
FLAGS="-g -Wall -Wextra -Werror -Wno-unused -Wno-unused-parameter -fsanitize=address"

gcc-7 -c iio.c -o iio.o $FLAGS
gcc-7 -c drawtriangle.c -o drawtriangle.o $FLAGS
gcc-7 -c drawtriangle.c -o drawtriangle.o $FLAGS
gcc-7 -c trimesh.c -o trimesh.o $FLAGS
gcc-7 -c rpc.c -o rpc.o $FLAGS
gcc-7 -c normals.c -o normals.o $FLAGS $(bash gdal-config --cflags)
gcc-7 -c get_utm_normal_shadow.c -o get_utm_normal_shadow.o $FLAGS $(bash gdal-config --cflags)
gcc-7 -c colorize_vertices_from_one_image.c -o colorize_vertices_from_one_image.o $FLAGS $(bash gdal-config --cflags)
gcc-7 -c write_coloured_ply.c -o write_coloured_ply.o $FLAGS $(bash gdal-config --cflags)
g++-7 -c geographiclib_wrapper.cpp -o geographiclib_wrapper.o $FLAGS $(bash gdal-config --cflags)

LIBS="-lz -ltiff -lpng -ljpeg -lm -lstdc++ -lGeographic -fsanitize=address"
gcc-7 iio.o drawtriangle.o trimesh.o rpc.o pickopt.o normals.o get_utm_normal_shadow.o geographiclib_wrapper.o -o ../bin/get_utm_normal_shadow $LIBS
gcc-7 iio.o trimesh.o rpc.o pickopt.o normals.o colorize_vertices_from_one_image.o geographiclib_wrapper.o -o ../bin/colorize_vertices_from_one_image $LIBS -framework GDAL
gcc-7 iio.o trimesh.o pickopt.o write_coloured_ply.o -o ../bin/write_coloured_ply $LIBS


