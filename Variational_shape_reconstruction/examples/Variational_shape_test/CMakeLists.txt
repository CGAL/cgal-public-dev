project (QEM)

cmake_minimum_required(VERSION 3.9)

# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

# Test if compiler has C++17
list(FIND CMAKE_CXX_COMPILE_FEATURES cxx_generalized_initializers has_cpp17)
if(has_cpp17 LESS 0)
  message(STATUS "NOTICE: These tests require a C++17 compiler and will not be compiled.")
  return()
endif()

# Use C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )

if( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

set(CGAL_DIR /usr/local/lib/CGAL)
set(Qt5_DIR /home/mathieu/lib_cgal/Qt5/5.12.12/gcc_64/lib/cmake/Qt5)
set(Qt5Charts_DIR /home/mathieu/lib_cgal/Qt5/5.12.12/gcc_64/lib/cmake/Qt5Charts)
set(OpenGL_GL_PREFERENCE GLVND) 

set(CMAKE_C_COMPILER "/usr/bin/clang")
set(CMAKE_CXX_COMPILER "/usr/bin/clang++")
set(OPENMP_LIBRARIES "/usr/local/Cellar/llvm/11.0.0/lib")
set(OPENMP_INCLUDES "/usr/local/Cellar/llvm/11.0.0/include")

set(OpenMP_CXX "${CMAKE_CXX_COMPILER}")
set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
set(OpenMP_CXX_LIB_NAMES "omp")
set(OpenMP_omp_LIBRARY omp)
# Find OpenMP
if(APPLE)
    set(CMAKE_C_COMPILER clang)
    set(CMAKE_CXX_COMPILER clang++)

    if(CMAKE_C_COMPILER_ID MATCHES "Clang\$")
        set(OpenMP_C "${CMAKE_C_COMPILER}")
        set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp")
        set(OpenMP_C_LIB_NAMES "omp")
        set(OpenMP_omp_LIBRARY omp)
    endif()

    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang\$")
        set(OpenMP_CXX "${CMAKE_CXX_COMPILER}")
        set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
        set(OpenMP_CXX_LIB_NAMES "omp")
        set(OpenMP_omp_LIBRARY omp)
    endif()

endif()

find_package(OpenMP REQUIRED)

find_package(CGAL COMPONENTS Qt5)
find_package(Qt5 QUIET COMPONENTS OpenGL)
find_package(Qt5Charts)
find_package(OpenGL)
find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)

find_package(SCIP QUIET)

set( SRCS main.cpp)

file(GLOB_RECURSE INCS "*.h")

if(WIN32)
	add_definitions( "/W3 /D_CRT_SECURE_NO_WARNINGS /wd4503 /nologo" )
endif(WIN32)

if(APPLE)
  add_compile_definitions(GL_SILENCE_DEPRECATION)
endif()


if( CGAL_FOUND AND Qt5_FOUND AND OPENGL_FOUND AND EIGEN3_FOUND AND OPENMP_FOUND AND SCIP_FOUND)
    include_directories(BEFORE . ./include ./ui)
    include(${CGAL_USE_FILE})
    include(${EIGEN3_USE_FILE})
    include_directories( BEFORE ${SCIP_INCLUDE_DIRS} )
    
    SET(MY_BUILD_PATH ${CMAKE_CURRENT_SOURCE_DIR}/build) #Variable build path
    SET(MY_SRC_PATH ${CMAKE_CURRENT_SOURCE_DIR})  #Variable src path
    INCLUDE_DIRECTORIES(${MY_BUILD_PATH})
    INCLUDE_DIRECTORIES(${MY_SRC_PATH})

    # UI files (Qt Designer files)
    qt5_wrap_ui( DT_UI_FILES ui/qem.ui)

    # The executable itself.
    add_executable( QEM ${SRCS} ${MOCS} ${DT_UI_FILES} ${INCS})
    qt5_use_modules(QEM OpenGL)

    # Link with Qt libraries
    target_link_libraries( QEM ${QT_LIBRARIES} )
    target_link_libraries( QEM Qt5::Charts )
    
    # Link with  OpenGL
    target_link_libraries( QEM ${OPENGL_LIBRARY})

    # Link with CGAL
    target_link_libraries( QEM ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
    
    # Link with SCIP
    target_link_libraries( QEM ${SCIP_LIBRARIES} )
    target_compile_definitions( QEM PUBLIC -DCGAL_USE_SCIP)

    # Link with openmp
    target_link_libraries( QEM OpenMP::OpenMP_CXX)

else()
  message(STATUS "NOTICE: This program requires CGAL, Qt5, OpenMP and Eigen. It will not be compiled.")

endif()
