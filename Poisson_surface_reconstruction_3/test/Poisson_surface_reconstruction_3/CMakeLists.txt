# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.15)
project( Implicit_surface_reconstruction_3_Tests )

# Find CGAL
find_package(CGAL REQUIRED)

# VisualC++ optimization for applications dealing with large data
if(MSVC)
  # Allow Windows 32bit applications to use up to 3GB of RAM
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /LARGEADDRESSAWARE")

  # Prints new compilation options
  message(
    STATUS
      "USING DEBUG CXXFLAGS   = '${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}'")
  message(
    STATUS
      "USING DEBUG EXEFLAGS   = '${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_DEBUG}'"
  )
  message(
    STATUS
      "USING RELEASE CXXFLAGS = '${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}'"
  )
  message(
    STATUS
      "USING RELEASE EXEFLAGS = '${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_RELEASE}'"
  )
endif()
  # Temporary debugging stuff
  find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)
  include(CGAL_Eigen3_support)
  if(TARGET CGAL::Eigen3_support)
    # Executables that require Eigen 3.1
    create_single_source_cgal_program("poisson_reconstruction_test.cpp")
    target_link_libraries(poisson_reconstruction_test PUBLIC CGAL::Eigen3_support)
    find_package(TBB)
    include(CGAL_TBB_support)
    if (TBB_FOUND)
      create_single_source_cgal_program( "poisson_and_parallel_mesh_3.cpp" )
      target_link_libraries(poisson_and_parallel_mesh_3 PUBLIC CGAL::Eigen3_support CGAL::TBB_support)
    else()
      message(STATUS "NOTICE: test with parallel Mesh_3 needs TBB and will not be compiled.")
    endif()

    find_package(Spectra REQUIRED)
    if(SPECTRA_FOUND)
      include_directories("${SPECTRA_INCLUDE_DIR}")
      find_package(Boost COMPONENTS system filesystem REQUIRED)
      add_definitions( "-BOOST_FILESYSTEM_LIBRARY" )
      add_definitions( "-BOOST_SYSTEM_LIBRARY" )
      add_definitions( "-DDEBUG_TRACE")
      list(APPEND CGAL_3RD_PARTY_LIBRARIES ${Boost_FILESYSTEM_LIBRARY})
      list(APPEND CGAL_3RD_PARTY_LIBRARIES ${Boost_SYSTEM_LIBRARY})
      create_single_source_cgal_program( "isr_test_mean_distance_mtp.cpp" )
      create_single_source_cgal_program( "isr_test_mean_distance_ptm.cpp" )
      create_single_source_cgal_program( "isr_test_hausdorff_ptm.cpp" )
      create_single_source_cgal_program( "isr_test_hausdorff_mtp.cpp" )
      create_single_source_cgal_program( "isr_test_mean_ang_dev.cpp" )
      create_single_source_cgal_program( "isr_test_check_genus.cpp" )
      create_single_source_cgal_program( "isr_test_hole_filling.cpp" )
      create_single_source_cgal_program( "isr_test_check_cc.cpp" )
      create_single_source_cgal_program( "isr_test_octree.cpp" )
    else()
      message(STATUS "NOTICE: Spectra was not found, tests won't be compiled.")
    endif()

  else()
    message(STATUS "NOTICE: Some of the executables in this directory need Eigen 3.1 (or greater) and will not be compiled.")
endif()
