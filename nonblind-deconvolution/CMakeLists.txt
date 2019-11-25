cmake_minimum_required(VERSION 2.8)
project(deblur)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -mtune=native -O3 -g")

include_directories(iio)

add_executable(deblur deblur.cpp image.cpp iio/iio.c)
set(LIBRARIES
    m
    fftw3_omp
    fftw3f_omp
    fftw3
    fftw3f
    pthread
    jpeg
    tiff
    png
)
target_link_libraries(deblur ${LIBRARIES})

