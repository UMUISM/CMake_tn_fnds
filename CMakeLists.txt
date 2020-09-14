cmake_minimum_required(VERSION 3.17)
project(CMake_tn_fnds)

set(CMAKE_CXX_STANDARD 11)

file(GLOB SRC src/*.cpp)
file(GLOB LIB lib/world/src/*.cpp)

include_directories(lib/world/src/)

add_library(WORLD ${LIB})

add_executable(CMake_tn_fnds ${SRC})
target_link_libraries(CMake_tn_fnds WORLD)