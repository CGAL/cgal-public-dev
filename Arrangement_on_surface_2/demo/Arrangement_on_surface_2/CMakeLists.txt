# This is the CMake script for compiling a CGAL application.

project( Arrangement_2 ) 

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
  cmake_policy(VERSION 2.8.4)
else()
  cmake_policy(VERSION 2.6)
endif()

find_package(CGAL COMPONENTS Core Qt4 )
include( ${CGAL_USE_FILE} )

find_package( Qt4 )

if ( CGAL_FOUND AND CGAL_Qt4_FOUND AND QT4_FOUND)
  if ( APPLE )
    add_definitions( -Wno-c++11-extensions )
  endif ( APPLE )

  # Uncomment this to enable verbose tracing of the sweep-line algorithm
  #add_definitions(-DCGAL_SL_VERBOSE )

  include( ${QT_USE_FILE} )
  include_directories( ./ )
  # Arrangement package includes
  include_directories( BEFORE ../../include )



  qt4_wrap_ui( arrangement_2_uis
    ArrangementDemoWindow.ui
    NewTabDialog.ui
    OverlayDialog.ui
    ArrangementDemoPropertiesDialog.ui
  )
  qt4_wrap_cpp( arrangement_2_mocs
    ArrangementDemoWindow.h
    ArrangementDemoPropertiesDialog.h
    ArrangementDemoTab.h
    Callback.h
    ColorItemEditor.h
    DeleteCurveModeItemEditor.h
    OverlayDialog.h
    PropertyValueDelegate.h
  )
  qt4_add_resources( arrangement_2_resources
    ArrangementDemoWindow.qrc
  )

  add_library( arr_demo_library SHARED
    ArrangementDemoGraphicsView.cpp
    ArrangementDemoPropertiesDialog.cpp
    ArrangementDemoTab.cpp
    ArrangementDemoWindow.cpp
    ArrangementGraphicsItem.cpp
    Callback.cpp
    ColorItemEditor.cpp
    CurveGraphicsItem.cpp
    DeleteCurveMode.cpp
    DeleteCurveModeItemEditor.cpp
    EnvelopeCallback.cpp
    FillFaceCallback.cpp
    GraphicsViewCurveInput.cpp
    NewTabDialog.cpp
    OverlayDialog.cpp
    PointsGraphicsItem.cpp
    PropertyValueDelegate.cpp
    SplitEdgeCallback.cpp
    Utils.cpp
    VerticalRayGraphicsItem.cpp
    VerticalRayShootCallback.cpp
    ${arrangement_2_mocs}
    ${arrangement_2_uis}
    ${arrangement_2_resources}
  )

  target_link_libraries( arr_demo_library
    ${CGAL_LIBRARIES}
    ${QT_LIBRARIES}
  )

  add_executable( arrangement_2
    arrangement_2.cpp

# provided by arr_demo_library
#    ArrangementDemoGraphicsView.cpp
#    ArrangementDemoPropertiesDialog.cpp
#    ArrangementDemoTab.cpp
#    ArrangementDemoWindow.cpp
#    ArrangementGraphicsItem.cpp
#    Callback.cpp
#    ColorItemEditor.cpp
#    DeleteCurveMode.cpp
#    DeleteCurveModeItemEditor.cpp
#    EnvelopeCallback.cpp
#    FillFaceCallback.cpp
#    GraphicsViewCurveInput.cpp
#    NewTabDialog.cpp
#    OverlayDialog.cpp
#    PointsGraphicsItem.cpp
#    PropertyValueDelegate.cpp
#    SplitEdgeCallback.cpp
#    Utils.cpp
#    VerticalRayGraphicsItem.cpp
#    VerticalRayShootCallback.cpp
#   ${arrangement_2_mocs}
#   ${arrangement_2_uis}
#   ${arrangement_2_resources}

#    GraphicsViewSegmentInput.cpp
  )

  target_link_libraries( arrangement_2
    arr_demo_library
    ${CGAL_LIBRARIES}
    ${CGAL_3RD_PARTY_LIBRARIES}
    ${QT_LIBRARIES}
  )

else()

  message(STATUS "NOTICE: This demo requires CGAL and Qt4, and will not be compiled.")

endif()

add_subdirectory( tests )
