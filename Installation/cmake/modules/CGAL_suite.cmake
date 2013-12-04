if(NOT CGAL_SUITE_FILE_INCLUDED)
  set(CGAL_SUITE_FILE_INCLUDED 1)

  include(CMakeParseArguments)
  include(CGAL_library)

  function(CGAL_are_depends_met DEPENDS)
    set(DEPENDS_MET TRUE)
    set(FAILED_DEPENDS)
    foreach(required_depend ${DEPENDS})
      if(${required_depend} MATCHES "CGAL")
        if(NOT TARGET ${required_depend})
          set(DEPENDS_MET FALSE)
          list(APPEND FAILED_DEPENDS ${required_depend})
        endif()
      elseif(NOT WITH_${required_depend})
        set(DEPENDS_MET FALSE)
        list(APPEND FAILED_DEPENDS ${required_depend})
      endif()
    endforeach()
    set(DEPENDS_MET ${DEPENDS_MET} PARENT_SCOPE)
    set(FAILED_DEPENDS ${FAILED_DEPENDS} PARENT_SCOPE)
  endfunction()
  
  function(CGAL_example)
    set(options) # no options
    set(oneValueArgs SOURCE TARGET)
    set(multiValueArgs REQUIRED_DEPENDENCIES OPTIONAL_DEPENDENCIES)

    cmake_parse_arguments(CGAL_example "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
    if("${CGAL_example_TARGET}" STREQUAL "")
      # strip the file ending
      string(REGEX REPLACE "\\..*$" "" CGAL_example_TARGET ${CGAL_example_SOURCE})
    endif()


    # Check if all required libraries are met.
    CGAL_are_depends_met("${CGAL_example_REQUIRED_DEPENDENCIES}")
    
    if(DEPENDS_MET)
      add_executable(${CGAL_example_TARGET} EXCLUDE_FROM_ALL ${CGAL_example_SOURCE})
      # Some examples use local include directories, support those.
      if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include)
        target_include_directories(${CGAL_example_TARGET} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
      endif()
      
      # the top-level dependency has the same name as the current project
      add_dependencies(${PROJECT_NAME} ${CGAL_example_TARGET})
      
      foreach(required_depend ${CGAL_example_REQUIRED_DEPENDENCIES})
        CGAL_use_library(${required_depend} ${CGAL_example_TARGET})
      endforeach()

      foreach(optional_depend ${CGAL_example_OPTIONAL_DEPENDENCIES})
        if(${lib} MATCHES "Boost") # We ignore WITH_ options for Boost.
          CGAL_use_library(optional_depend ${source_clean})
        elseif(${lib} MATCHES "CGAL" AND TARGET ${lib}) # We have a target for this internal dependency.
          CGAL_use_library(optional_depend ${source_clean})
        elseif(WITH_${lib})
          CGAL_use_library(optional_depend ${source_clean})
        endif()
      endforeach()
    else()
      message(STATUS "The example ${CGAL_example_SOURCE} requires ${FAILED_DEPENDS} and is not being build.")
    endif()

  endfunction()

  macro(CGAL_example_suite)
    set(options) # no options
    set(oneValueArgs NAME)
    set(multiValueArgs REQUIRED_DEPENDENCIES OPTIONAL_DEPENDENCIES SOURCES)
    cmake_parse_arguments(CGAL_example_suite "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
    project(${CGAL_example_suite_NAME} CXX)

    # Check if all required libraries are met. We don't leave this to
    # CGAL_example, so we can emit a better diagnostic.
    CGAL_are_depends_met("${CGAL_example_suite_REQUIRED_DEPENDENCIES}")

    if(DEPENDS_MET)
      add_custom_target(${PROJECT_NAME})

      if(TARGET examples) # Don't do this if we are a in a stand-alone build.
        add_dependencies(examples ${PROJECT_NAME})
      else()
        set_property(TARGET ${PROJECT_NAME} PROPERTY EXCLUDE_FROM_ALL OFF)
      endif()

      foreach(source ${CGAL_example_suite_SOURCES})
        CGAL_example(
          REQUIRED_DEPENDENCIES ${CGAL_example_suite_REQUIRED_DEPENDENCIES}
          OPTIONAL_DEPENDENCIES ${CGAL_example_suite_OPTIONAL_DEPENDENCIES}
          SOURCE ${source})
      endforeach()
    else()
      message(STATUS "The examples ${CGAL_example_suite_SOURCES} require ${FAILED_DEPENDS} and are not being build.")
    endif()

  endmacro()
endif()