project (CGAL_Symbolic)
   
# Minimal version of CMake:
cmake_minimum_required(VERSION 2.6.2)

# Tested version:
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

if ( NOT CGAL_FOUND )   
    find_package(CGAL QUIET COMPONENTS Core)
endif()

include(${CGAL_USE_FILE})

if ( NOT MPFI_FOUND )
# optionally look for support for gmp number types
find_package( MPFI QUIET) 
if ( MPFI_FOUND )
 include( ${MPFI_USE_FILE} )
endif( MPFI_FOUND )

endif()

# try to find CUDA

INCLUDE (CheckCXXSourceRuns)
# INCLUDE (CheckCSourceCompiles)

# TODO EBEB add WITH_GPU_DIR
# TODO EBEB environment variables do not work, you should use
# find_path
# TODO EBEB add CUDA_LIB_DIR
# set(CUDA_TOOLKIT_ROOT_DIR "$ENV{HOME}/cuda3.0/cuda")

# set ENV variables:
# export CUDA_BIN_PATH=$HOME/cuda3.0/cuda/bin
# export CUDA_LIB_PATH=$HOME/cuda3.0/cuda/lib

find_package(CUDA QUIET REQUIRED)

if ( CUDA_FOUND )

  set(CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE OFF)
  set(CUDA_VERBOSE_BUILD ON)
  set(CUDA_64_BIT_DEVICE_CODE OFF)

  set(CMAKE_REQUIRED_INCLUDES ${CUDA_INCLUDE_DIRS})
  set(CMAKE_REQUIRED_LIBRARIES ${CUDA_LIBRARIES})
  FILE(READ test_CUDA.C src)

  CHECK_CXX_SOURCE_RUNS("${src}" CUDA_CHECK)

  if ( CUDA_CHECK_COMPILED )

    message(STATUS  "test_CUDA compiled successfully")
    MATH(EXPR CUDA_MAJOR '${CUDA_CHECK_EXITCODE}/10')
    MATH(EXPR CUDA_MINOR '${CUDA_CHECK_EXITCODE}%10')

    if (CUDA_CHECK_EXITCODE LESS 13)

      message(STATUS "\nDevice does not support CUDA or version is too small\n")

    else()

      set(GPU_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR})
# "/KM/home-0/emeliyan/work/CGAL/branches/features/Symbolic-mpi/Symbolic/src/CGALSymbolic")

      if(CUDA_CHECK_EXITCODE LESS 20)
          cache_set(CGAL_Symbolic_3RD_PARTY_LIBRARIES   ${CUDA_LIBRARIES}
              "${GPU_LIB_DIR}/libgpu_symbolic_sm_13.so")
      else()
          cache_set(CGAL_Symbolic_3RD_PARTY_LIBRARIES   ${CUDA_LIBRARIES}
              "${GPU_LIB_DIR}/libgpu_symbolic_sm_20.so")
      endif()

      cache_set(CGAL_Symbolic_3RD_PARTY_DEFINITIONS  "-DCGAL_USE_GPU")
      cache_set(CGAL_Symbolic_USE_GPU "ON")

      message(STATUS "\nCUDA version found: ${CUDA_CHECK_EXITCODE}\n")

    endif()

  else( CUDA_CHECK_COMPILED )

    message(STATUS  "test_CUDA failed: ${OUTPUT}")

  endif( CUDA_CHECK_COMPILED )

endif( CUDA_FOUND )

# try to find optional NTL

if ( NOT NTL_FOUND )
  find_package( NTL )
endif( NOT NTL_FOUND )

if ( NTL_FOUND )

# TODO   get_dependency_version(NTL NTL)

  message( STATUS "NTL include:     ${NTL_INCLUDE_DIR}" )
  message( STATUS "NTL libraries:   ${NTL_LIBRARIES}" )

  cache_set(CGAL_Symbolic_3RD_PARTY_INCLUDE_DIRS
        ${CGAL_Symbolic_3RD_PARTY_INCLUDE_DIRS} ${NTL_INCLUDE_DIR} )
  cache_set(CGAL_Symbolic_3RD_PARTY_LIBRARIES
        ${CGAL_Symbolic_3RD_PARTY_LIBRARIES} ${NTL_LIBRARIES} )
  cache_set(CGAL_Symbolic_3RD_PARTY_DEFINITIONS
        ${CGAL_Symbolic_3RD_PARTY_DEFINITIONS} "-DCGAL_USE_NTL")
  cache_set(CGAL_Symbolic_USE_NTL "ON")

endif ( NTL_FOUND )

include_directories( ${CGAL_3RD_PARTY_INCLUDE_DIRS} ${CGAL_Symbolic_3RD_PARTY_INCLUDE_DIRS} )
  
build_cgal_library( CGAL_Symbolic CGALSymbolic "")

# add_library(CGAL_Symbolic base.cpp)

add_dependencies( CGAL_Symbolic CGAL )
  
add_definitions( ${CGAL_3RD_PARTY_DEFINITIONS} ${CGAL_Symbolic_3RD_PARTY_DEFINITIONS} )
  
set_target_properties(CGAL_Symbolic PROPERTIES VERSION "${CGAL_SOVERSION}" SOVERSION "${CGAL_SONAME_VERSION}")
  
target_link_libraries( CGAL_Symbolic ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} ${CGAL_Symbolic_3RD_PARTY_LIBRARIES} )
  
set( CGAL_Symbolic_BASENAME CGAL_Symbolic )

if (NOT CGAL_AUTO_LINK_ENABLED )
    
  set_target_properties( CGAL_Symbolic PROPERTIES OUTPUT_NAME "${CGAL_Symbolic_BASENAME}" )
    
  get_target_property( CGAL_Symbolic_LIBRARY_FULL_PATH CGAL_Symbolic LOCATION  )
    
  set(CGAL_Symbolic_LIBRARY "${CGAL_Symbolic_LIBRARY_FULL_PATH}" CACHE STRING "The CGAL_Symbolic library" FORCE)
    
  mark_as_advanced(CGAL_Symbolic_LIBRARY)
    
  target_link_libraries( CGAL_Symbolic ${CGAL_LIBRARY} )
    
else()
    
  set(CGAL_Symbolic_LIBRARY)
    
  set(CGAL_Symbolic_RELEASE_MANGLED_NAME        "${CGAL_Symbolic_BASENAME}-${CGAL_TOOLSET}-mt" ) 
  set(CGAL_Symbolic_DEBUG_MANGLED_NAME          "${CGAL_Symbolic_BASENAME}-${CGAL_TOOLSET}-mt-gd" ) 
  set(CGAL_Symbolic_MINSIZEREL_MANGLED_NAME     "${CGAL_Symbolic_BASENAME}-${CGAL_TOOLSET}-mt-o" ) 
  set(CGAL_Symbolic_RELWITHDEBINFO_MANGLED_NAME "${CGAL_Symbolic_BASENAME}-${CGAL_TOOLSET}-mt-g" ) 
    
  set_target_properties( CGAL_Symbolic PROPERTIES DEBUG_OUTPUT_NAME          "${CGAL_Symbolic_DEBUG_MANGLED_NAME}" 
    RELEASE_OUTPUT_NAME        "${CGAL_Symbolic_RELEASE_MANGLED_NAME}" 
    MINSIZEREL_OUTPUT_NAME     "${CGAL_Symbolic_MINSIZEREL_MANGLED_NAME}" 
    RELWITHDEBINFO_OUTPUT_NAME "${CGAL_Symbolic_RELWITHDEBINFO_MANGLED_NAME}" 
    )
    
  if ( HAS_CFG_INTDIR )                      
    add_custom_command(TARGET CGAL_Symbolic POST_BUILD COMMAND copy /Y \"$(TargetPath)\" \"$(SolutionDir)lib\" )
    add_custom_command(TARGET CGAL_Symbolic POST_BUILD COMMAND if exist \"$(TargetDir)$(TargetName).pdb\" copy /Y \"$(TargetDir)$(TargetName).pdb\" \"$(SolutionDir)lib\" )
  endif()
    
endif()
  
if ( CGAL_INSTALL_LIB_DIR )
  install(TARGETS CGAL_Symbolic DESTINATION ${CGAL_INSTALL_LIB_DIR} )
endif()
    
