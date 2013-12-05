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
      
    else(NOT ${lib_pos} EQUAL -1)
      if(DEFINED CGAL_WITH_${lib})
        set(use_it ${CGAL_WITH_${lib}})
      else()
        set(use_it ${${vlib}_FOUND})
      endif()

      if(NOT use_it)
        return()
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
    elseif(${lib_pos} EQUAL -1)
      message(FATAL_ERROR "CGAL_use_library called with ${lib} which is not a supported library.")
    endif()
  endfunction()
endif()