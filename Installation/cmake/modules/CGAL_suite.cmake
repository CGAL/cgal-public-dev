if(NOT CGAL_SUITE_FILE_INCLUDED)
  set(CGAL_SUITE_FILE_INCLUDED 1)

  include(CMakeParseArguments)
  include(CGAL_library)

  function(CGAL_are_depends_met DEPENDS)
    set(DEPENDS_MET TRUE)
    set(FAILED_DEPENDS)
    foreach(required_depend ${DEPENDS})
      if(${required_depend} MATCHES "CGAL")
        if(NOT BUILD_${required_depend})
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
    set(oneValueArgs SOURCE)
    set(multiValueArgs REQUIRED_DEPENDENCIES OPTIONAL_DEPENDENCIES)

    cmake_parse_arguments(CGAL_example "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    # Check if all required libraries are met.
    CGAL_are_depends_met(${CGAL_example_REQUIRED_DEPENDENCIES})

    if(DEPENDS_MET)
      # strip the file ending
      string(REGEX REPLACE "\\..*$" "" source_clean ${CGAL_example_SOURCE})
      add_executable(${source_clean} EXCLUDE_FROM_ALL ${source})
      foreach(required_depend ${CGAL_example_REQUIRED_DEPENDENCIES})
        CGAL_use_library(${required_depend} ${source_clean})
      endforeach()

      foreach(optional_depend ${CGAL_example_OPTIONAL_DEPENDENCIES})
        # CGAL_use_library(required_depend ${source_clean})
      endforeach()
    else()
      message(STATUS "The examples ${CGAL_example_suite_SOURCES} require
        ${FAILED_DEPENDS} and are not being build.")
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
    CGAL_are_depends_met(${CGAL_example_suite_REQUIRED_DEPENDENCIES})

    if(DEPENDS_MET)
      foreach(source ${CGAL_example_suite_SOURCES})
        CGAL_example(
          REQUIRED_DEPENDENCIES ${CGAL_example_suite_REQUIRED_DEPENDENCIES}
          OPTIONAL_DEPENDENCIES ${CGAL_example_suite_OPTIONAL_DEPENDENCIES}
          SOURCE ${source})
      endforeach()
    else()
      message(STATUS "The examples ${CGAL_example_suite_SOURCES} require
        ${FAILED_DEPENDS} and are not being build.")
    endif()

  endmacro()
endif()