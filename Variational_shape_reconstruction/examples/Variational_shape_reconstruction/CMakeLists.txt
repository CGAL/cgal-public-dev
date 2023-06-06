cmake_minimum_required(VERSION 3.1...3.23)
project(Variational_shape_reconstruction_Examples)

find_package(CGAL REQUIRED)

create_single_source_cgal_program("variational_shape_reconstruction_example.cpp")