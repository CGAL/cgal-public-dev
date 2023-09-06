cmake_minimum_required( VERSION 3.1...3.23 )
project( gsoc2023 )

set( CMAKE_CXX_EXTENSIONS OFF )

include( CheckIPOSupported )
check_ipo_supported( RESULT ipo_supported )
if( ipo_supported )
  set( CMAKE_INTERPROCEDURAL_OPTIMIZATION True )
endif()

add_subdirectory(examples/Continuous_collision_detection_3)
