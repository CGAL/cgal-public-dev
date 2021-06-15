cmake_minimum_required(VERSION 3.13)
project(ray-bbox-intersections-array_of_structs)
set(CMAKE_CXX_STANDARD 17)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -Wall -march=native")

add_executable(array_of_structs
        array_of_structs/main.cpp)

add_executable(struct_of_arrays
        struct_of_arrays/main.cpp)

add_executable(array_of_batched_structs
        array_of_batched_structs/main.cpp)
