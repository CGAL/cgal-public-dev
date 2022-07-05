cmake_minimum_required(VERSION 3.1...3.23)
project(OctreeMarchingCubes)

find_package(CGAL REQUIRED QUIET OPTIONAL_COMPONENTS Core)
if (CGAL_FOUND)

    create_single_source_cgal_program("main.cpp")

    find_package(Eigen3 3.1.91) #(requires 3.2.0 or greater)
    include(CGAL_Eigen_support)

else ()
    message(WARNING
            "This program requires the CGAL library, and will not be compiled.")
endif ()
