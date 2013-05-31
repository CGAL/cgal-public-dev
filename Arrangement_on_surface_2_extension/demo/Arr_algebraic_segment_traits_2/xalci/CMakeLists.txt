# Created by the script cgal_create_cmake_script_with_options
# This is the CMake script for compiling a set of CGAL applications.

project( xalci )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

# set(CMAKE_MODULE_PATH "$ENV{HOME}/work/CGAL/branches/unsorted-branches/CMAKE/Installation/cmake/modules")
# /home/emeliyan/work/CGAL/branches/unsorted-branches/CMAKE/Installation/cmake/modules:

set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )
 
if ( COMMAND cmake_policy )

  cmake_policy( SET CMP0003 NEW )  

endif()

# NOTE NOTE: it is important that Qt3 linker flags appear *after* CUDA
# libraries otherwise runtime crash
# CGAL and its components
find_package( CGAL REQUIRED COMPONENTS Core Symbolic Qt3 )

if ( NOT CGAL_FOUND )
  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  
endif()

# include helper file
include( ${CGAL_USE_FILE} )

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

message(STATUS "PATH: ${CMAKE_MODULE_PATH} ")

# Boost and its components
# find_package( Boost REQUIRED )
# 
# if ( NOT Boost_FOUND )
# 
#   message(STATUS "This project requires the Boost library, and will not be compiled.")
# 
#   return()  
# 
# endif()

# link_libraries( ${Boost_THREAD_LIBRARY} )

# include for local directory
include_directories( BEFORE
$ENV{HOME}/work/CGAL/branches/unsorted-branches/asm/Arrangement_on_surface_2_extension/include $ENV{HOME}/work/CGAL/branches/experimental-packages/Algebraic_kernel_d/include include  )

# include for local package

# include of additional packages
# include_directories( thread/include/ )

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


# Creating entries for target: xalci
# ############################

if ( CGAL_Qt3_FOUND AND QT3_FOUND )

qt3_automoc(  misc.cpp xalci.cpp )

  # Make sure the compiler can find generated .moc files
  include_directories( BEFORE ${CMAKE_CURRENT_BINARY_DIR} )
 
  include_directories( ${QT3_INCLUDE_DIR} )

endif()

add_executable( xalci  misc.cpp xalci.cpp )

# target_link_libraries(xalci ${CGAL_3RD_PARTY_LIBRARIES} ${CGAL_LIBRARIES})

add_to_cached_list( CGAL_EXECUTABLE_TARGETS xalci )

message(STATUS "------- CGAL LIBS: ${CGAL_LIBRARIES}")
message(STATUS "------- CGAL_3rd_party: ${CGAL_3RD_PARTY_LIBRARIES}")

link_libraries( ${Boost_THREAD_LIBRARY} "$ENV{HOME}/thirdparty_libs/ntl32/lib/libntl.a")

# Link the executable to CGAL and third-party libraries
# target_link_libraries(xalci  ${QT3_LIBRARIES} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

# target_link_libraries(xalci ${CGAL_3RD_PARTY_LIBRARIES} ${CGAL_LIBRARIES})


