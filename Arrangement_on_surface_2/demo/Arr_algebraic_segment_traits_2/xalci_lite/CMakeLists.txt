
project( xalci_lite )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )
 
if ( COMMAND cmake_policy )

  cmake_policy( SET CMP0003 NEW )  

endif()

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m32")

find_package( CGAL REQUIRED COMPONENTS Core Qt3 )

if ( NOT CGAL_FOUND )
  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  
endif()

# include helper file
include( ${CGAL_USE_FILE} )

# find_package( MPFI )
# if ( MPFI_FOUND )
#  include( ${MPFI_USE_FILE} )
# else()
#   message(STATUS "This project requires the MPFI library, and will not be compiled.")
#   return()
# endif( MPFI_FOUND ) 
# 
# find_package( RS )
# if ( RS_FOUND )
#  include( ${RS_USE_FILE} )
# else()
#   message(STATUS "This project requires the RS library, and will not be compiled.")
#   return()
# endif( RS_FOUND )   
# 

# Boost and its c   omponents
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()

endif()

link_libraries( ${Boost_THREAD_LIBRARY} )

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

  # Make sure the compiler can find generated .moc files
  include_directories( ${CMAKE_CURRENT_BINARY_DIR} include )
 
  include_directories( ${QT3_INCLUDE_DIR} )

  qt3_automoc(  misc.cpp xalci.cpp )

endif()

add_executable( xalci_lite  misc.cpp xalci.cpp )

include_directories(include $ENV{HOME}/work/CGAL/branches/unsorted-branches/Numerical_algebraic_kernel_d/include $ENV{HOME}/work/CGAL/branches/experimental-packages/Algebraic_kernel_2/include
$ENV{HOME}/work/CGAL/branches/experimental-packages/Algebraic_kernel_d/include
$ENV{HOME}/work/CGAL/branches/unsorted-branches/Arcavoid/include $ENV{HOME}/work/CGAL/branches/features/Arrangement_on_surface_2/include)


add_to_cached_list( CGAL_EXECUTABLE_TARGETS xalci_lite )

link_libraries( ${Boost_THREAD_LIBRARY})


