if(NOT CGAL_LIBRARY_FILE_INCLUDED)
  set(CGAL_LIBRARY_FILE_INCLUDED 1)

  include(CMakeParseArguments)
  include(CMakeDependentOption)
  include(CGAL_use_library)

  # Define a build target for a library.
  #
  # Example:
  #
  # CGAL_define_library(
  #   NAME Foobar # required
  #   BUILD_DEFAULT OFF # should this library be built by default?
  #   REQUIRED_DEPENDENCIES CGAL GMP # optional, can be external or internal dependencies
  #   OPTIONAL_DEPENDENCIES Eigen3   # optional, can be external or internal dependencies
  #   SOURCES foobar0.cpp foobar1.cpp foobar2.cpp # optional, if empty sources are collected by globbing
  # )
  macro(CGAL_define_library)
    set(options) # no options
    set(oneValueArgs NAME BUILD_DEFAULT)
    set(multiValueArgs REQUIRED_DEPENDENCIES OPTIONAL_DEPENDENCIES SOURCES)
    cmake_parse_arguments(CGAL_define_library "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
    
    message(STATUS "Configuring lib${CGAL_define_library_NAME}")
    list(APPEND CGAL_CONFIGURED_LIBRARIES ${CGAL_define_library_NAME})
    set(CGAL_CONFIGURED_LIBRARIES ${CGAL_CONFIGURED_LIBRARIES} PARENT_SCOPE)

    option(BUILD_${CGAL_define_library_NAME} "Build the ${CGAL_define_library_NAME} library"
      ${CGAL_define_library_BUILD_DEFAULT})

    if(BUILD_${CGAL_define_library_NAME})
      # define the target
      list(LENGTH CGAL_define_library_SOURCES CGAL_define_library_source_length)
      if(CGAL_define_library_source_length EQUAL 0)
        foreach (package ${CGAL_CONFIGURED_PACKAGES})
          file(GLOB CGAL_LIBRARY_SOURCE_FILES_TMP ${package}/src/${CGAL_define_library_NAME}/*.cpp)
          list(APPEND CGAL_LIBRARY_SOURCE_FILES ${CGAL_LIBRARY_SOURCE_FILES_TMP})
        endforeach()
        add_library(${CGAL_define_library_NAME} ${CGAL_LIBRARY_SOURCE_FILES})
        # collect_cgal_library(${CGAL_define_library_NAME} "")
      else()
        # rather broken
        add_library(${CGAL_define_library_NAME} ${CGAL_define_library_SOURCES})
      endif()

      # set the CGAL include directories
      target_include_directories(${CGAL_define_library_NAME} PUBLIC ${CGAL_INCLUDE_DIRS})

      foreach(required_depend ${CGAL_define_library_REQUIRED_DEPENDENCIES})
        set(BUILD_${CGAL_define_library_NAME}_WITH_${required_depend} TRUE)
        CGAL_use_library(${required_depend} ${CGAL_define_library_NAME})
      endforeach()

      foreach(opt_depend ${CGAL_define_library_OPTIONAL_DEPENDENCIES})
        # Only show this option if we are actually building this library
        # and if the dependency is enabled.
        CMAKE_DEPENDENT_OPTION(BUILD_${CGAL_define_library_NAME}_WITH_${opt_depend}
          "Build library ${CGAL_define_library_NAME} with ${opt_depend} support" ON
          "BUILD_${CGAL_define_library_NAME};WITH_${opt_depend}" OFF)
        if(BUILD_${CGAL_define_library_NAME}_WITH_${opt_depend})
          CGAL_use_library(${opt_depend} ${CGAL_define_library_NAME})
        endif()
      endforeach()

      # export the library target
      export(TARGETS ${CGAL_define_library_NAME}
        APPEND FILE ${CGAL_EXPORT_FILE})
    endif()
  endmacro()
  
endif()