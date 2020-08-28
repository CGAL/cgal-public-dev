# This is the CMake script for compiling this folder.

cmake_minimum_required(VERSION 3.1...3.15)
project( Poisson_surface_reconstruction_3_Examples )


# Find CGAL
find_package(CGAL QUIET)

if ( CGAL_FOUND )

  # VisualC++ optimization for applications dealing with large data
  if (MSVC)
    # Allow Windows 32bit applications to use up to 3GB of RAM
    SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /LARGEADDRESSAWARE")

    # Print new compilation options
    message( STATUS "USING DEBUG CXXFLAGS   = '${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}'" )
    message( STATUS "USING DEBUG EXEFLAGS   = '${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_DEBUG}'" )
    message( STATUS "USING RELEASE CXXFLAGS = '${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}'" )
    message( STATUS "USING RELEASE EXEFLAGS = '${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_RELEASE}'" )
  endif()

  # Find Eigen3 (requires 3.1.0 or greater)
  find_package(Eigen3 3.1.0)
  if (EIGEN3_FOUND)
    # Executables that require Eigen 3
    create_single_source_cgal_program( "poisson_reconstruction_example.cpp" )
    CGAL_target_use_Eigen(poisson_reconstruction_example)
    create_single_source_cgal_program( "poisson_reconstruction.cpp" )
    CGAL_target_use_Eigen(poisson_reconstruction)
    create_single_source_cgal_program( "poisson_reconstruction_function.cpp" )
    CGAL_target_use_Eigen(poisson_reconstruction_function)
    create_single_source_cgal_program( "tutorial_example.cpp" )
    CGAL_target_use_Eigen(tutorial_example)
  else()
    message(STATUS "NOTICE: The examples need Eigen 3.1 (or greater) will not be compiled.")
  endif()


else()

    message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")

endif()



# NEEDS FOR ADDING IN-MESHING
set(source_dir "/home/felix/code/source")
include_directories(${source_dir})

# Options
# =====================================================================================================================
option(WITH_POISSON_REF   "Compile the reference Screened Poisson Reconstruction" OFF)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_property(GLOBAL PROPERTY PREDEFINED_TARGETS_FOLDER "cmake")

# Dependencies & utility
# =====================================================================================================================
if(WITH_POISSON_REF)
  find_package(OpenMP REQUIRED)
endif()

function(target_compile_to_bin target)
  set_target_properties(${target} PROPERTIES
          VS_DEBUGGER_WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/bin
          RUNTIME_OUTPUT_DIRECTORY      $<1:${CMAKE_CURRENT_LIST_DIR}/bin>
          ARCHIVE_OUTPUT_DIRECTORY      $<1:${CMAKE_CURRENT_LIST_DIR}/bin>
          LIBRARY_OUTPUT_DIRECTORY      $<1:${CMAKE_CURRENT_LIST_DIR}/bin>)
endfunction()
function(adjust_exe_target target)
  source_group(${source_dir} FILES ${ARGN})
  target_compile_to_bin(${target})
  target_link_libraries(${target} PUBLIC core)
endfunction()

# Poisson (reference)
# =====================================================================================================================
if(WITH_POISSON_REF)
  file(GLOB REF_SOURCES ${source_dir}/poisson_recon/*.h ${source_dir}/poisson_recon/*.inl)
  source_group("poisson_recon" FILES ${REF_SOURCES} ${source_dir}/poisson_recon/PoissonRecon.cpp ${source_dir}/poisson_recon/PlyFile.cpp)
  add_executable(poisson_recon_ref ${REF_SOURCES} ${source_dir}/poisson_recon/PoissonRecon.cpp ${source_dir}/poisson_recon/PlyFile.cpp)
  target_link_libraries(poisson_recon_ref PUBLIC OpenMP::OpenMP_CXX)
  if(MSVC)
    target_compile_options(poisson_recon_ref PUBLIC /arch:AVX2)
  endif()
  target_compile_to_bin(poisson_recon_ref)
endif()

# Poisson: custom implementation
# =====================================================================================================================
file(GLOB_RECURSE SOURCES ${source_dir}/poisson_surface/*.h ${source_dir}/poisson_surface/*.cpp ${source_dir}/poisson_surface/*.inl)
source_group(TREE ${source_dir} FILES ${SOURCES})
add_library(poisson ${SOURCES})
CGAL_target_use_Eigen(poisson)
target_include_directories(poisson PUBLIC ${source_dir})
target_compile_features(poisson PUBLIC cxx_std_17)
set_target_properties(poisson PROPERTIES FOLDER "core")
if(MSVC)
  target_compile_definitions(poisson PUBLIC _SILENCE_ALL_CXX17_DEPRECATION_WARNINGS)
endif()

# In-meshing: core
# =====================================================================================================================
file(GLOB_RECURSE SOURCES ${source_dir}/core/*.h ${source_dir}/core/*.cpp ${source_dir}/core/*.inl)
set(SOURCES ${SOURCES} ${source_dir}/poisson_recon/Bridge.cpp ${source_dir}/poisson_recon/Bridge.h ${source_dir}/base.h ${source_dir}/debugbreak.h)
source_group(TREE ${source_dir} FILES ${SOURCES})
message("SOURCES : ")
message(${SOURCES})
add_library(core ${SOURCES})
target_link_libraries(core PUBLIC poisson)
set_target_properties(core PROPERTIES FOLDER "core")

# In-meshing: command line
# =====================================================================================================================
create_single_source_cgal_program( "reconstruction_with_inmeshing.cpp" )
CGAL_target_use_Eigen(reconstruction_with_inmeshing)
#target_link_libraries(reconstruction_with_inmeshing PUBLIC core)

adjust_exe_target(reconstruction_with_inmeshing reconstruction_with_inmeshing.cpp)
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT main)

# NEEDS FOR ADDING IN-MESHING - END