cmake_minimum_required(VERSION 3.27)
project(Mat)

set(CMAKE_CXX_STANDARD 20)

include_directories(include)


find_package(Boost REQUIRED)

add_subdirectory(examples)