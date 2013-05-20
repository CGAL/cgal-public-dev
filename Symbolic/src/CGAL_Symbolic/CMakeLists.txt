message("Configuring libCGAL_Symbolic")

# External libraries needed

if ( NOT MPFI_FOUND )
  # optionally look for support for gmp number types
  find_package( MPFI QUIET)
  if ( MPFI_FOUND )
    include( ${MPFI_USE_FILE} )
   endif( MPFI_FOUND )
endif()

# TODO what if MPFI is not found?

# try to find CUDA

INCLUDE (CheckCXXSourceRuns)

set(CUDA_TOOLKIT_ROOT_DIR "$ENV{CUDA_TOOLKIT_ROOT_DIR}")

find_package(CUDA QUIET)

if  ( CUDA_FOUND )

# TODO TODO TODO: compile lib with debug support

#   set(CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE OFF)
#   set(CUDA_VERBOSE_BUILD ON)
#   set(CUDA_64_BIT_DEVICE_CODE OFF)

  set(CMAKE_REQUIRED_INCLUDES ${CUDA_INCLUDE_DIRS})
  set(CMAKE_REQUIRED_LIBRARIES ${CUDA_LIBRARIES})
  FILE(READ test_CUDA.C src)

  CHECK_CXX_SOURCE_RUNS("${src}" CUDA_CHECK)

  if ( CUDA_CHECK_COMPILED )

    message(STATUS  "test_CUDA compiled successfully: ${CUDA_CHECK_EXITCODE}")
    MATH(EXPR CUDA_MAJOR '${CUDA_CHECK_EXITCODE}/10')
    MATH(EXPR CUDA_MINOR '${CUDA_CHECK_EXITCODE}%10')

    message("exit code: ${CUDA_CHECK_EXITCODE}")

    if (CUDA_CHECK_EXITCODE LESS 13)

      message("")
      message("   WARNING: Device does not support CUDA or version is too small")
      message("")

    else()

      # TODO EBEB Find GPU_LIB_DIR by FindGPUSYM (which would imply GPU_LIB_DIR to be set as env-variable)
      set(GPU_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR})

      set(libEXT ${CMAKE_SHARED_LIBRARY_SUFFIX}) 

      if(CUDA_CHECK_EXITCODE LESS 20)
        if(CMAKE_SIZEOF_VOID_P EQUAL 8)
            cache_set(CGAL_Symbolic_3RD_PARTY_LIBRARIES  
                 ${GPU_LIB_DIR}/libgpu_symbolic_sm_13_64${libEXT}  ${CUDA_LIBRARIES})
        else()
            cache_set(CGAL_Symbolic_3RD_PARTY_LIBRARIES  
                 ${GPU_LIB_DIR}/libgpu_symbolic_sm_13_32${libEXT}  ${CUDA_LIBRARIES})  
        endif()
      else()
        if(CMAKE_SIZEOF_VOID_P EQUAL 8)
            cache_set(CGAL_Symbolic_3RD_PARTY_LIBRARIES  
                 ${GPU_LIB_DIR}/libgpu_symbolic_sm_20_64${libEXT}  ${CUDA_LIBRARIES})
        else()
            cache_set(CGAL_Symbolic_3RD_PARTY_LIBRARIES  
                 ${GPU_LIB_DIR}/libgpu_symbolic_sm_20_32${libEXT}  ${CUDA_LIBRARIES})  
        endif()
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

if(COMMAND add_config_flag)
  set( CGAL_HAS_SYMBOLIC TRUE ) 
  add_config_flag( CGAL_HAS_SYMBOLIC )
endif()

use_essential_libs()

include_directories( ${CGAL_3RD_PARTY_INCLUDE_DIRS} ${CGAL_Symbolic_3RD_PARTY_INCLUDE_DIRS} )
  
add_definitions( ${CGAL_3RD_PARTY_DEFINITIONS} ${CGAL_Symbolic_3RD_PARTY_DEFINITIONS} )
 
link_directories    ( ${CGAL_LIBRARIES_DIR} ${CGAL_3RD_PARTY_LIBRARIES_DIRS} )
 
collect_cgal_library( CGAL_Symbolic "")

add_dependencies( CGAL_Symbolic CGAL )
  
set_target_properties(CGAL_Symbolic PROPERTIES VERSION "${CGAL_SOVERSION}" SOVERSION "${CGAL_SONAME_VERSION}")
  
target_link_libraries( CGAL_Symbolic ${CGAL_LIBRARY} ${CGAL_Core_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} ${CGAL_Symbolic_3RD_PARTY_LIBRARIES} )

# TODO what about the following lines? Other libs' CMakeLists.txt do
# not contain these lined

  
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
    
message("libCGAL_Symbolic is configured")
