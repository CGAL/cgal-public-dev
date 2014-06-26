add_executable( Test Test.cpp )
add_executable( TestRayPL TestRayPL.cpp )
add_executable( TestAlgebraic TestAlgebraic.cpp )
include_directories( BEFORE ${CMAKE_CURRENT_BINARY_DIR} )
include_directories( BEFORE ../../ )
include_directories( BEFORE ../ )
message( STATUS ${CMAKE_CURRENT_BINARY_DIR} )

set( bezier_example_srcs
  BezierExampleWindow.cpp
)
qt4_wrap_ui( bezier_example_uis
  BezierExampleWindow.ui
)
qt4_wrap_cpp( bezier_example_mocs
  BezierExampleWindow.h
)
qt4_add_resources( bezier_example_resources
  ../ArrangementDemoWindow.qrc
)
add_executable( bezier_example bezier_example.cpp
  ${bezier_example_srcs}
  ${bezier_example_uis}
  ${bezier_example_mocs}
  ${bezier_example_resources}
  ../Utils.cpp
  ../ArrangementGraphicsItem.cpp
  ../ArrangementDemoGraphicsView.cpp
  ../PointSetGraphicsItem.cpp
)
target_link_libraries( bezier_example
    ${CGAL_LIBRARIES}
)

set( bezier_example_2_srcs
  BezierTabWindow.cpp
  ../Callback.cpp
  ../ArrangementDemoTab.cpp
  ../ArrangementGraphicsItem.cpp
  ../VerticalRayShootCallback.cpp
  ../EnvelopeCallback.cpp
  ../SplitEdgeCallback.cpp
  ../FillFaceCallback.cpp
  ../GraphicsViewCurveInput.cpp
  ../VerticalRayGraphicsItem.cpp
  ../PointsGraphicsItem.cpp
  ../Utils.cpp
  ../ArrangementDemoGraphicsView.cpp

)
qt4_wrap_ui( bezier_example_2_uis
  BezierTabWindow.ui
)
qt4_wrap_cpp( bezier_example_2_mocs
  BezierTabWindow.h
  ../Callback.h
  ../ArrangementDemoTab.h
)
add_executable( bezier_example_2 bezier_example_2.cpp
  ${bezier_example_2_srcs}
  ${bezier_example_2_uis}
  ${bezier_example_2_mocs}
)

target_link_libraries( bezier_example_2
  ${CGAL_LIBRARIES}
  ${QT_LIBRARIES}
)

add_executable( bezier_envelope bezier_envelope.cpp )
