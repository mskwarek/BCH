cmake_minimum_required(VERSION 3.5)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

file(GLOB SRC_FILES *.cpp)
add_library(BCH_codes ${SRC_FILES} Decoder.cpp Decoder.h Encoder.cpp Encoder.h GfField.h Polynomial.h)