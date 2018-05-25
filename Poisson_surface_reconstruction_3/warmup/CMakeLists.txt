project(mesh)

cmake_minimum_required(VERSION 2.8.11)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

if(POLICY CMP0043)
    cmake_policy( SET CMP0043  OLD )
endif()

# Include this package's headers first
include_directories( BEFORE ./ ./include ../../include )

# Add specific Find.cmake modules
set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules )



find_package(CGAL COMPONENTS Qt5)
include( ${CGAL_USE_FILE} )

find_package(Qt5 QUIET COMPONENTS Xml OpenGL)
find_package(OpenGL)

if( Qt5_FOUND )
    find_package(QGLViewer)
endif()

add_definitions(-DBOOST_NO_CXX11_TEMPLATE_ALIASES)

#set( HDRS c2t3.h options.h viewer.h scene.h window.h function_ellipsoid.h function_sphere.h function_torus.h types.h)
set( HDRS c2t3.h options.h viewer.h scene.h window.h function_ellipsoid.h function_sphere.h function_torus.h types.h)

set( SRCS viewer.cpp main.cpp scene.cpp window.cpp )

set( MOCS moc_window.cxx
          moc_options.cxx
          moc_viewer.cxx )

if( CGAL_Qt5_FOUND AND Qt5_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND )
    include_directories( ${QGLVIEWER_INCLUDE_DIR} )

    # UI files (Qt Designer files)
    Qt5_wrap_ui( DT_UI_FILES mesh.ui options.ui)

    include(AddFileDependencies)

    # qrc files (resources files, that contain icons, at least)
    Qt5_add_resources( DT_RESOURCE_FILES mesh.qrc )

    # use the Qt MOC preprocessor on classes that derives from QObject
    Qt5_generate_moc( window.h moc_window.cxx )
    Qt5_generate_moc( viewer.h moc_viewer.cxx )
    Qt5_generate_moc( options.h moc_options.cxx )

    # The executable itself.
    add_executable( mesh ${SRCS} ${MOCS} ${DT_UI_FILES} ${DT_RESOURCE_FILES} ${HDRS})

    qt5_use_modules(mesh OpenGL Xml)

    # Link with Qt libraries
    target_link_libraries( mesh ${QT_LIBRARIES} )

    # Link with CGAL
    target_link_libraries( mesh ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

    # Link with libQGLViewer, OpenGL
    target_link_libraries( mesh ${QGLVIEWER_LIBRARIES} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} )

    add_to_cached_list( CGAL_EXECUTABLE_TARGETS mesh )

else()

    set(OT_MISSING_DEPS "")

    if(NOT OPENGL_FOUND)
        set(OT_MISSING_DEPS "OpenGL, ${OT_MISSING_DEPS}")
    endif()

    if(NOT Qt5_FOUND)
        set(OT_MISSING_DEPS "Qt5, ${OT_MISSING_DEPS}")
    endif()

    if(NOT QGLVIEWER_FOUND)
        set(OT_MISSING_DEPS "QGLViewer, ${OT_MISSING_DEPS}")
    endif()

    if(NOT CGAL_FOUND)
        set(OT_MISSING_DEPS "the CGAL library, ${OT_MISSING_DEPS}")
    endif()

    if(NOT CGAL_Qt5_FOUND)
        set(OT_MISSING_DEPS "the CGAL Qt5 library, ${OT_MISSING_DEPS}")
    endif()

    message(STATUS "NOTICE: This demo requires ${OT_MISSING_DEPS} and will not be compiled.")

endif()
