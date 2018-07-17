cmake_minimum_required(VERSION 3.6)
project(GC)

# Are we using gcc?
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # GCC on MacOs needs this option to use the clang assembler
    if (APPLE)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wa,-q")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-q")
    endif ()
    # Optimize to the current CPU and enable warnings
    #	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native -std=c99 -O3")
    #	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -O3")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -g")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
endif ()

# Default build is in Release mode
IF(NOT CMAKE_BUILD_TYPE AND NOT MSVC)
    SET(CMAKE_BUILD_TYPE "Release")
ENDIF(NOT CMAKE_BUILD_TYPE AND NOT MSVC)
SET(CMAKE_BUILD_TYPE "Debug")

SET (LIBS "")

FIND_PACKAGE(TIFF)
IF(TIFF_FOUND)
    INCLUDE_DIRECTORIES(${TIFF_INCLUDE_DIR})
    #  LINK_DIRECTORIES(${TIFF_LIBRARY})
    SET(LIBS ${LIBS} ${TIFF_LIBRARIES})
ENDIF()

FIND_PACKAGE(JPEG)
if(JPEG_FOUND)
    INCLUDE_DIRECTORIES(${JPEG_INCLUDE_DIR})
    # LINK_DIRECTORIES(${JPEG_LIBRARY})
    SET(LIBS ${LIBS} ${JPEG_LIBRARIES})
ENDIF()

FIND_PACKAGE(PNG)
IF(PNG_FOUND)
    INCLUDE_DIRECTORIES(${PNG_INCLUDE_DIR})
    #  LINK_DIRECTORIES(${PNG_LIBRARY})
    SET(LIBS ${LIBS} ${PNG_LIBRARIES})
ENDIF()


# DETECT LINUX
IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    # Linux specific code flags
    SET(C_FLAGS "${C_FLAGS} -DI_CAN_HAS_LINUX -D_GNU_SOURCE")
ENDIF()


# IIO
SET(C_FLAGS "${C_FLAGS} -O3 -std=c99 -funroll-loops -Wno-unused -DNDEBUG")

#FFTW
find_path (FFTW_INCLUDES fftw3.h)
find_library (FFTW_LIBRARIES NAMES fftw3)
include_directories(GC ${FFTW_INCLUDES})

# Source files
set(SOURCE_FILES iio.c iio.h interpolation.h interpolation.cpp gaussFilter.h gradient.h gradient.cpp matops.h matops.cpp gradcorr.h gradcorr.cpp window.h window.cpp satreg.h satreg.cpp)

add_executable(registerAerial ${SOURCE_FILES} registerAerial.cpp)
target_link_libraries(registerAerial ${FFTW_LIBRARIES})
target_link_libraries(registerAerial ${LIBS})

add_executable(gc ${SOURCE_FILES} gc.cpp)
target_link_libraries(gc ${FFTW_LIBRARIES})
target_link_libraries(gc ${LIBS})

add_executable(registerAll ${SOURCE_FILES} registerAll.cpp)
target_link_libraries(registerAll ${FFTW_LIBRARIES})
target_link_libraries(registerAll ${LIBS})

