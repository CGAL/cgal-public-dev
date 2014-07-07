add_executable( Test Test.cpp )
add_executable( TestRayPL TestRayPL.cpp )
add_executable( TestAlgebraic TestAlgebraic.cpp )
include_directories( BEFORE ${CMAKE_CURRENT_BINARY_DIR} )
include_directories( BEFORE ../../ )
include_directories( BEFORE ../ )
message( STATUS ${CMAKE_CURRENT_BINARY_DIR} )

add_executable( test_polynomial_io test_polynomial_io.cpp )

#####
# Test code for scaffolding the Bezier graphics item.
#####
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

#####
# Bezier-only demo.
#####
set( bezier_example_2_srcs
  BezierTabWindow.cpp
#  ../Callback.cpp
#  ../ArrangementDemoTab.cpp
#  ../ArrangementGraphicsItem.cpp
#  ../VerticalRayShootCallback.cpp
#  ../EnvelopeCallback.cpp
#  ../SplitEdgeCallback.cpp
#  ../FillFaceCallback.cpp
#  ../GraphicsViewCurveInput.cpp
#  ../VerticalRayGraphicsItem.cpp
#  ../PointsGraphicsItem.cpp
#  ../Utils.cpp
#  ../ArrangementDemoGraphicsView.cpp
)
qt4_wrap_ui( bezier_example_2_uis
  BezierTabWindow.ui
)
qt4_wrap_cpp( bezier_example_2_mocs
  BezierTabWindow.h
#  ../Callback.h
#  ../ArrangementDemoTab.h
)
add_executable( bezier_example_2 bezier_example_2.cpp
  ${bezier_example_2_srcs}
  ${bezier_example_2_uis}
  ${bezier_example_2_mocs}
)

target_link_libraries( bezier_example_2
  arr_demo_library
  ${CGAL_LIBRARIES}
  ${QT_LIBRARIES}
)

#####
# Small example for computing the envelope of a Bezier arrangement.
#####
add_executable( bezier_envelope bezier_envelope.cpp )

#####
# Small demo compiled for one particular arrangement geometry trait.
#####
add_library( arr_demo_traits_library SHARED
  BezierDemoTraits.cpp
  AlgebraicDemoTraits.cpp
)

set( single_tab_example_srcs
  SingleTabWindow.cpp
)
qt4_wrap_ui( single_tab_example_uis
  SingleTabWindow.ui
)
qt4_wrap_cpp( single_tab_example_mocs
  SingleTabWindow.h
)
add_executable( single_tab_example single_tab_example.cpp
  ${single_tab_example_srcs}
  ${single_tab_example_uis}
  ${single_tab_example_mocs}
)
target_link_libraries( single_tab_example
  arr_demo_library
  arr_demo_traits_library
  ${CGAL_LIBRARIES}
  ${QT_LIBRARIES}
)
