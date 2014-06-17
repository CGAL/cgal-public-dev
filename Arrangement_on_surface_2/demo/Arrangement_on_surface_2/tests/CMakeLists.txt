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
add_executable( bezier_example bezier_example.cpp
  ${bezier_example_srcs}
  ${bezier_example_uis}
  ${bezier_example_mocs}
  ../Utils.cpp
  ../ArrangementGraphicsItem.cpp
  ../ArrangementDemoGraphicsView.cpp
)
target_link_libraries( bezier_example
    ${CGAL_LIBRARIES}
)
