cmake_minimum_required(VERSION 3.13)
project(ray-bbox-intersections-array_of_structs)
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -Wall -march=native")

find_package(CGAL)

add_executable(main
        main.cpp)

target_link_libraries(main CGAL::CGAL)