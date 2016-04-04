# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Envelope_voronoi_2_examples )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
 
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()
 
find_package(CGAL QUIET COMPONENTS Core Qt3 )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )
  
  find_package(Qt3-patched QUIET )
  # FindQt3-patched.cmake is FindQt3.cmake patched by CGAL developers, so
  # that it can be used together with FindQt4: all its variables are prefixed
  # by "QT3_" instead of "QT_".
  
  if(CGAL_Qt3_FOUND AND QT3_FOUND)

    include( Qt3Macros-patched )
    qt3_automoc(  triangle_area_distance_voronoi.cpp )
    qt3_automoc(  mobius_diagram.cpp )
    
    # Link with Boost.ProgramOptions
    find_package(Boost QUIET COMPONENTS program_options)
    if(Boost_PROGRAM_OPTIONS_FOUND)
      if( CGAL_AUTO_LINK_ENABLED )
	message( STATUS "Boost.ProgramOptions library: found" )
      else()
	message( STATUS "Boost.ProgramOptions library: ${Boost_PROGRAM_OPTIONS_LIBRARY}" )
      endif()
      add_definitions( "-DCGAL_USE_BOOST_PROGRAM_OPTIONS" )
      if ( NOT CGAL_AUTO_LINK_ENABLED )
	link_libraries( ${Boost_PROGRAM_OPTIONS_LIBRARY} )
      endif()

      # Make sure the compiler can find generated .moc files
      include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR})
      
      
      include_directories(../../include)
      include_directories(../../../../../trunk/Envelope_3/include)
      include_directories(../../../../../trunk/Polynomial/include)
      include_directories(../../../../../trunk/Number_types/include)
      include_directories(../../../../../trunk/Algebraic_foundations/include)
      include_directories(../../../../../trunk/Modular_arithmetic/include)
      include_directories(../../../Arithmetic_kernel/include)
      include_directories(../../../Curved_kernel_via_analysis_2/include)
      include_directories(../../../Voronoi_diagram_of_lines_3/include)
      include_directories(../../../Algebraic_kernel_d/include)
      include_directories(../../../mpfr/include)

      include_directories( ${QT3_INCLUDE_DIR} )


      add_executable  (triangle_area_distance_voronoi  triangle_area_distance_voronoi.cpp )
      add_executable  (mobius_diagram  mobius_diagram.cpp )
      add_executable  (power_diagram_on_sphere power_diagram_on_sphere.cpp )
      add_executable  (voronoi_of_points_on_sphere voronoi_of_points_on_sphere.cpp )
      add_executable  (anisotropic_diagram anisotropic_diagram.cpp )
      
      add_to_cached_list( CGAL_EXECUTABLE_TARGETS triangle_area_distance_voronoi )
      add_to_cached_list( CGAL_EXECUTABLE_TARGETS mobius_diagram )
      add_to_cached_list( CGAL_EXECUTABLE_TARGETS power_diagram_on_sphere )
      add_to_cached_list( CGAL_EXECUTABLE_TARGETS voronoi_of_points_on_sphere )
      add_to_cached_list( CGAL_EXECUTABLE_TARGETS anisotropic_diagram )
        
      # Link the executable to CGAL and third-party libraries
      link_libraries(${QT3_LIBRARIES} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

  else(Boost_PROGRAM_OPTIONS_FOUND)
    
    message(STATUS "NOTICE: This example requires Boost Program Options, and will not be compiled.")
    
  endif()
else()
  
    message(STATUS "NOTICE: This example requires Qt3 and the CGAL Qt3 library, and will not be compiled.")
  
  endif()
  
else()
  
    message(STATUS "NOTICE: This example requires the CGAL library, and will not be compiled.")
  
endif()

