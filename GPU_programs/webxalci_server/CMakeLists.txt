
##############################################################################
# There are two ways to add include directories to the NVCC command
# line:

cmake_minimum_required(VERSION 2.8)

project(webxalci_server)

#set(CUDA_SDK_ROOT_DIR "$ENV{HOME}/cuda/sdk/C")
# set(CUDA_INCLUDE_DIRS "$ENV{HOME}/cuda/sdk/C/common/inc")

# set(CMAKE_MODULE_PATH "$ENV{CMAKE_MODULE_PATH}")
set(CMAKE_BUILD_TYPE Release)
# set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m32 -O2")

find_package(CGAL QUIET COMPONENTS MPFI RS3 Core Symbolic Qt3)

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


find_package( RS3 )
if ( RS3_FOUND )
 include( ${RS3_USE_FILE} )
else()
  message(STATUS "This project requires the RS library, and will not be compiled.")
  return()
endif( RS3_FOUND )

# find_package( CGAL REQUIRED COMPONENTS Symbolic)

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

# qt3_automoc(  misc.cpp xalci.cpp )

  # Make sure the compiler can find generated .moc files
  include_directories( BEFORE ${CMAKE_CURRENT_BINARY_DIR} )
 
  include_directories( ${QT3_INCLUDE_DIR} )

endif()

## NOTE NOTE NOTE it is important to place
## 'features/Symbolic-mpi/Algebraic_kernel_d' *AFTER*
## 'Numerical_algebraic_kernel_d' since the latter one has its own copy of RS
## which, otherwise, generates duplicate symbols if there are several C sources

## NOTE NOTE: if there are problems with MPFI: edit CGAL_UseMPFI.cmake in
## Installation directory

include_directories(BEFORE $ENV{HOME}/work/CGAL/branches/experimental-packages/Algebraic_kernel_d/include $ENV{HOME}/work/CGAL/branches/unsorted-branches/Arcavoid/include
$ENV{HOME}/work/CGAL/branches/features/Arrangement_on_surface_2/include $ENV{HOME}/thirdparty_libs/ntl32/include .)

# NOTE NOTE NOTE: server needs to be compiled with -lcuda (driver API)
# also make sure that your version of CUDA SDK matches that of installed
# driver API (which otherwise can cause serious problems)

# NOTE NOTE NOTE: CUDA-related libraries should be linked *before* OpenGL and
# probably Qt stuff, otherwise program crashes at runtime !!!
# As of now, there is no way to guarantee this..

# set(USE_32_BIT_MODULI_SET true)
# 
# if(USE_32_BIT_MODULI_SET)
# 
# set(GPU_LIBS "$ENV{SYMBLIB}/src/CGALSymbolic/libgpu_symbolic_sm_20.so"
#     "$ENV{HOME}/cuda3.0/cuda/lib/libcudart.so" "-lcuda" )
# 
# else()
# 
# set(GPU_LIBS "$ENV{SYMBLIB}/src/CGALSymbolic/libgpu_symbolic_sm_13.so"
#      "$ENV{HOME}/cuda3.0/cuda/lib/libcudart.so" "-lcuda")
# 
# endif()

# link_libraries( ${GPU_LIBS} )

 message(STATUS "cgal 3rdparty libs: ${CGAL_3RD_PARTY_LIBRARIES}")

 message(STATUS "cgal libs : ${CGAL_LIBRARIES}") 

add_executable( webxalci_server rasterizer.C server.C shm_buffer.C helper_functions.C )

target_link_libraries(webxalci_server  "$ENV{HOME}/thirdparty_libs/ntl32/lib/libntl.a" -lssl ${CGAL_3RD_PARTY_LIBRARIES} ${CGAL_LIBRARIES} -lCGAL_Core)

# "$ENV{HOME}/thirdparty_libs/ntl32/lib/libntl.a"
# target_link_libraries( ${Boost_THREAD_LIBRARY}  -lssl)

# Add a special target to clean nvcc generated files.
# CUDA_BUILD_CLEAN_TARGET()

