if(NOT CGAL_LIBRARY_FILE_INCLUDED)
  set(CGAL_LIBRARY_FILE_INCLUDED 1)

  include(CMakeParseArguments)

  # Define a build target for a library.
  #
  # Example:
  #
  # CGAL_define_library(
  #   NAME Foobar # required
  #   REQUIRED_DEPENDENCIES CGAL GMP # optional, can be external or internal dependencies
  #   OPTIONAL_DEPENDENCIES Eigen3   # optional, can be external or internal dependencies
  #   SOURCES foobar0.cpp foobar1.cpp foobar2.cpp # optional, if empty sources are collected by globbing
  # )
  macro(CGAL_define_library)
    set(options) # no options
    set(oneValueArgs NAME)
    set(multiValueArgs REQUIRED_DEPENDENCIES OPTIONAL_DEPENDENCIES SOURCES)
    cmake_parse_arguments(CGAL_define_library "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
    message("Configuring lib${CGAL_define_library_NAME}")

    option(BUILD_${CGAL_define_library_NAME} "Build the ${CGAL_define_library_name} library")

    foreach(required_depend ${CGAL_define_library_REQUIRED_DEPENDENCIES})
      set(BUILD_${CGAL_define_library_NAME}_WITH_${required_depend} TRUE)
      CGAL_use_library(${required_depend})
    endforeach()

    foreach(opt_depend ${CGAL_define_library_OPTIONAL_DEPENDENCIES})
      option(BUILD_${CGAL_define_library_NAME}_WITH_${opt_depend}
        "Build library ${CGAL_define_library_NAME} with ${opt_depend} support" FALSE)
      if(BUILD_${CGAL_define_library_NAME}_WITH_${opt_depend})
        CGAL_use_library(${opt_depend})
      endif()
    endforeach()

    list(LENGTH CGAL_define_library_SOURCES CGAL_define_library_source_length)
    if(CGAL_define_library_source_length EQUAL 0)
      collect_cgal_library(${CGAL_define_library_NAME} "")
    else()
      # rather broken
      add_library(${CGAL_define_library_NAME}
        ${CGAL_define_library_SOURCES})
    endif()

    # target_link_libraries(${CGAL_define_library_NAME})
    message("lib${CGAL_define_library_NAME} is configured")
  endmacro()

endif()