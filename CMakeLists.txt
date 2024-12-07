cmake_minimum_required(VERSION 3.25)
project(Mat)

set(CMAKE_CXX_STANDARD 23)

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE "-march=native -Ofast -g")

include_directories(include)

add_subdirectory(examples)
add_subdirectory(test)
add_subdirectory(benchmark)