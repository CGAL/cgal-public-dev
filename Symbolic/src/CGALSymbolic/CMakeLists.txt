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
  find_package(CGAL REQUIRED)
endif()

include(${CGAL_USE_FILE})

# TODO NTL optional
# TODO GPU optional

include_directories( ${CGAL_3RD_PARTY_INCLUDE_DIRS} ${CGAL_Symbolic_3RD_PARTY_INCLUDE_DIRS} )
  
build_cgal_library( CGAL_Symbolic CGALSymbolic "")
  
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
    
