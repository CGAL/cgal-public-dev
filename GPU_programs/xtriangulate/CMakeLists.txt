##############################################################################
# There are two ways to add include directories to the NVCC command
# line:

cmake_minimum_required(VERSION 2.8)

project(xtriangulate)

set(CUMP_COMPILE_LIBRARY 0)

#set(CUDA_SDK_ROOT_DIR "$ENV{HOME}/cuda/sdk/C")
# set(CUDA_INCLUDE_DIRS "$ENV{HOME}/cuda/sdk/C/common/inc")

# set(CMAKE_MODULE_PATH "$ENV{CMAKE_MODULE_PATH}")
# set(CMAKE_BUILD_TYPE Release)
set(CMAKE_BUILD_TYPE Debug)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m32")
set(CMAKE_CXX_FLAGS "-O2 -g -m32")

find_package(CGAL QUIET COMPONENTS MPFI RS Core Symbolic Qt3)

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()

# include helper file
include( ${CGAL_USE_FILE} )

find_package( Boost REQUIRED COMPONENTS thread )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()

endif()

find_package( MPFI )
if ( MPFI_FOUND )
 include( ${MPFI_USE_FILE} )
else()
  message(STATUS "This project requires the MPFI library, and will not be compiled.")
  return()
endif( MPFI_FOUND )


find_package( RS )
if ( RS_FOUND )
 include( ${RS_USE_FILE} )
else()
  message(STATUS "This project requires the RS library, and will not be compiled.")
  return()
endif( RS_FOUND )

find_package( CGAL REQUIRED COMPONENTS Symbolic)

# Qt3
# FindQt3-patched.cmake is FindQt3.cmake patched by CGAL developers, so
# that it can be used together with FindQt4: all its variables are prefixed
# by "QT3_" instead of "QT_".
find_package(Qt3-patched QUIET )

if ( NOT QT3_FOUND )

  message(STATUS "This project requires the Qt3 library, and will not be compiled.")
  return()  

endif()

if ( CGAL_Qt3_FOUND )
  
  include( Qt3Macros-patched )

endif()

if ( CGAL_Qt3_FOUND AND QT3_FOUND )

  # Make sure the compiler can find generated .moc files
  include_directories( ${CMAKE_CURRENT_BINARY_DIR} include )
 
  include_directories( ${QT3_INCLUDE_DIR} )

    qt3_automoc(  mainwnd.cpp ogl_view.cpp )

endif()

## NOTE NOTE NOTE it is important to place
## 'features/Symbolic-mpi/Algebraic_kernel_d' *AFTER*
## 'Numerical_algebraic_kernel_d' since the latter one has its own copy of RS
## which, otherwise, generates duplicate symbols if there are several C sources

# NOTE NOTE NOTE: 'include' must be specified for 'include_directories' otherwise qt3_automoc cannot find the headers to generate moc files

# $ENV{HOME}/work/CGAL/branches/experimental-packages/Algebraic_kernel_2/include

include_directories(BEFORE include $ENV{HOME}/work/CGAL/branches/experimental-packages/Algebraic_kernel_d/include $ENV{HOME}/work/CGAL/branches/unsorted-branches/Arcavoid/include $ENV{HOME}/work/CGAL/branches/features/Arrangement_on_surface_2/include
$ENV{HOME}/work/CGAL/branches/experimental-packages/Curved_kernel_via_analysis_2_plus_1/include $ENV{HOME}/work/CGAL/branches/experimental-packages/Algebraic_kernel_3/include $ENV{HOME}/thirdparty_libs/ntl32/include )


# NOTE NOTE NOTE: server needs to be compiled with -lcuda (driver API)
# also make sure that your version of CUDA SDK matches that of installed
# driver API (which otherwise can cause serious problems)

# NOTE NOTE NOTE: CUDA-related libraries should be linked *before* OpenGL and
# probably Qt stuff, otherwise program crashes at runtime !!!
# As of now, there is no way to guarantee this..

add_executable( xtriangulate ogl_view.cpp mainwnd.cpp zy_skeletonizer.cpp zz_surf_analysis.cpp xtriangulate.cpp ogl_exts.cpp math.cpp)

target_link_libraries(xtriangulate "-lglut" "$ENV{HOME}/thirdparty_libs/ntl32/lib/libntl.a")

# target_link_libraries(xtriangulate   )


