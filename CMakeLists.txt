cmake_minimum_required(VERSION 3.5)
project(BCH)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp)

include_directories(${PROJECT_SOURCE_DIR}/BCH_codes)

add_subdirectory(BCH_codes)
add_subdirectory(bch_tests)

link_libraries(BCH_codes)
add_executable(BCH ${SOURCE_FILES})