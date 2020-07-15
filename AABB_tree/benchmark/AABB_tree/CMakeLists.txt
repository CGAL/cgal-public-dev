# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.14)
project( AABB_tree_benchmark )

find_package(CGAL  QUIET )
FIND_PACKAGE(embree 3.0 QUIET)
find_package(TBB REQUIRED)

include(CGAL_TBB_support)

include_directories(include)

if ( CGAL_FOUND AND embree_FOUND )

    ADD_EXECUTABLE(AabbTest AabbTest.cpp include/RaysGenerate.h)
    TARGET_LINK_LIBRARIES(AabbTest CGAL )

    ADD_EXECUTABLE(EmbreeTest EmbreeTest.cpp include/RaysGenerate.h)
    TARGET_LINK_LIBRARIES(EmbreeTest CGAL embree)

if ( TBB_FOUND )
    ADD_EXECUTABLE(EmbreeParallel EmbreeParallel.cpp include/RaysGenerate.h)
    TARGET_LINK_LIBRARIES(EmbreeParallel PUBLIC CGAL embree)

    ADD_EXECUTABLE(AabbParallel AabbParallel.cpp include/RaysGenerate.h)
    TARGET_LINK_LIBRARIES(AabbParallel PUBLIC CGAL  CGAL::TBB_support)
endif()

    ADD_EXECUTABLE(UserGeomSM UserGeomSM.cpp include/RaysGenerate.h)
    TARGET_LINK_LIBRARIES(UserGeomSM CGAL embree)

else()

    message(STATUS "This program requires the CGAL and the Embree library, and will not be compiled.")

endif()
