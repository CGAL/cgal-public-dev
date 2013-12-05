if(NOT CGAL_USE_LIBRARY_FILE_INCLUDED)
  set(CGAL_USE_LIBRARY_FILE_INCLUDED 1)

  function(CGAL_use_library lib target_name)
    set(optional FALSE)
    set(vlib ${CGAL_${lib}_PREFIX})
    if(TARGET ${lib}) 
      # This is true for all CGAL libraries and for imported targets.
      # Takes care of all transitive dependencies, include
      # directories, and definitions.
      target_link_libraries(${target_name} ${lib}) 
    elseif(${lib} MATCHES "Qt")
      
    else()
      # confirm we support it
      list(FIND CGAL_EXTERNAL_LIBRARIES ${lib} lib_pos)
      if(${lib_pos} EQUAL -1)
        message(FATAL_ERROR "${lib} is not a supported external library.")
      endif()

      if(NOT ${vlib}_FOUND AND NOT optional)
        message(FATAL_ERROR "${target_name} requires ${lib} but it could not be found.")
      endif()
    endif()

    
    if(NOT WITH_${lib})
      message(FATAL_ERROR "${target_name} requires ${lib}, but WITH_${lib} is \"OFF\".")
    endif()
    
    if(NOT ${vlib}_FOUND)
      message(FATAL_ERROR "Trying to use ${lib} which could not be found.")
    endif()

    # if none of the above set use_define, set it here
    if(NOT DEFINED use_define)
      set(use_define ${vlib})
    endif()
    
    if(${lib} MATCHES "Boost_") # a boost component library has to link with singular LIBRARY
      target_link_libraries(${target_name} ${${vlib}_LIBRARY})
    else()
      target_link_libraries(${target_name} ${${vlib}_LIBRARIES})
    endif()
    
    if(${CMAKE_VERSION} VERSION_LESS 2.8.12)
      target_include_directories(${target_name} PUBLIC "${${vlib}_INCLUDE_DIR}")
    else()
      target_include_directories(${target_name} SYSTEM PUBLIC "${${vlib}_INCLUDE_DIR}")
    endif()
    
    target_compile_definitions(${target_name} PUBLIC "${${vlib}_DEFINITIONS}" PUBLIC "-DCGAL_USE_${use_define}")
  endfunction()

endif()