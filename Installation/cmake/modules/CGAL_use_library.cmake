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

      if(${lib} MATCHES "Boost_") # a boost component library has to
                                  # link with singular LIBRARY and
                                  # cannot use normal vlib_include_dir
        target_link_libraries(${target_name} ${${vlib}_LIBRARY})
        if(${CMAKE_VERSION} VERSION_LESS 2.8.12)
          target_include_directories(${target_name} PUBLIC "${Boost_INCLUDE_DIR}")
        else()
          target_include_directories(${target_name} SYSTEM PUBLIC "${Boost_INCLUDE_DIR}")
        endif()
      else()
        target_link_libraries(${target_name} ${${vlib}_LIBRARIES})
        if(${CMAKE_VERSION} VERSION_LESS 2.8.12)
          target_include_directories(${target_name} PUBLIC "${${vlib}_INCLUDE_DIR}")
        else()
          target_include_directories(${target_name} SYSTEM PUBLIC "${${vlib}_INCLUDE_DIR}")
        endif()
      endif()
      
      # Bug: this adds e.g. CGAL_USE_Boost_SYSTEM while we use
      # CGAL_USE_BOOST_SYSTEM in code. What to fix?
      target_compile_definitions(${target_name} PUBLIC "${${vlib}_DEFINITIONS}" PUBLIC "-DCGAL_USE_${vlib}")
    elseif(${lib_pos} EQUAL -1)
      message(FATAL_ERROR "CGAL_use_library called with ${lib} which is not a supported library.")
    endif()
  endfunction()
endif()