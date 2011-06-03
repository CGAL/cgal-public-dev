# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Kernel_d_luis_example )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
 
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()
 
find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE ../../include)
  #include_directories (BEFORE /home/luis/SVN/linbox)
  include_directories (BEFORE /home/luis/src/linbox-1.2.0)
  include_directories (BEFORE /home/luis/src/fflas-ffpack-1.4.1)
  #include_directories (BEFORE /home/luis/src/linbox-1.1.7)
  #include_directories (BEFORE /opt/linbox-1.1.7-gmp4/include)
  #include_directories (BEFORE /opt/givaro-3.3.5/include)
  #include_directories (BEFORE /opt/givaro-3.3.5-gmp4/include)
  include_directories (BEFORE /opt/givaro-trunk/include)
  include_directories (BEFORE /opt/TOPCOM-0.16.2/include)

  add_definitions(${CMAKE_CXX_FLAGS} "-DLinBoxSrcOnly")
  #add_definitions(${CMAKE_CXX_FLAGS} "-DDISABLE_COMMENTATOR")
  #link_directories(/opt/givaro-3.3.5/lib)
  #link_directories(/opt/givaro-3.3.5-gmp4/lib)
  link_directories(/opt/givaro-trunk/lib)
  #link_libraries("blas")
  link_libraries("lapack")
  link_libraries("givaro")

  create_single_source_cgal_program( "linbox_traits.cpp" )
  create_single_source_cgal_program( "linbox_mpq_class_traits.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

