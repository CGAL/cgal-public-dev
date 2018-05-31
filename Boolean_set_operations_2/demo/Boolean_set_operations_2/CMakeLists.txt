cmake_minimum_required(VERSION 2.8)

project(Boolean_set_operations_2)


cmake_minimum_required(VERSION 2.8.11)
if(POLICY CMP0043)
  cmake_policy(SET CMP0043 OLD)
endif()


# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

# Require c++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(CGAL COMPONENTS Qt5 Core )

include(${CGAL_USE_FILE})

set(QT_USE_QTMAIN TRUE)

find_package(Qt5 QUIET COMPONENTS Widgets Script Svg OpenGL)

include_directories (BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/HeaderFiles )

include_directories (BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/../../include )
include_directories (BEFORE ${CMAKE_CURRENT_BINARY_DIR})
#include_directories (BEFORE ${CMAKE_CURRENT_BINARY_DIR}/include)


qt5_add_resources(UI_RESOURCES
    UiFiles/Boolean_operations_2.qrc)

qt5_wrap_ui(UI_HEADERS
    UiFiles/Boolean_set_operations_2.ui
	UiFiles/ChoiseDialogPolygon.ui)

add_executable(Boolean_set_operations_2
	Boolean_set_operations_2.cpp
    #MainWindow.cpp

    #HeaderFiles/MainWindow.h
    HeaderFiles/QT5/Gps_segement_traits_2_apurva.h
	HeaderFiles/Typedefs.h
	HeaderFiles/QT5/BoundaryPiecesGraphicsItem.h
	HeaderFiles/QT5/CircularPolygons.h
	HeaderFiles/QT5/GraphicsViewCircularPolygonInput.h
	HeaderFiles/QT5/PiecewiseBoundaryGraphicsItem.h
	HeaderFiles/QT5/PiecewiseGraphicsItemBase.h
	HeaderFiles/QT5/PiecewiseRegionGraphicsItem.h
	HeaderFiles/QT5/PiecewiseSetGraphicsItem.h
	HeaderFiles/QT5/GraphicsViewLinearPolygonInput.h
	HeaderFiles/QT5/LinearPolygons.h
    ${UI_HEADERS}
    ${UI_RESOURCES})

# Link with Qt libraries
target_link_libraries(Boolean_set_operations_2 ${QT_LIBRARIES})

# Link with CGAL
target_link_libraries(Boolean_set_operations_2 ${CGAL_LIBRARIES}  ${CGAL_3RD_PARTY_LIBRARIES})


qt5_use_modules(Boolean_set_operations_2
    Widgets
    Script
    Svg
	OpenGL)

#set_target_properties(Boolean_set_operations_2 PROPERTIES LINK_FLAGS "/SUBSYSTEM:WINDOWS /ENTRY:mainCRTStartup")

#	HeaderFiles/QT5/GpsCircleSegmentGraphicsItem.h
# 	HeaderFiles/GuiDrawPolygon.h
#	SourceFiles/ChoiseDialogPolygon.cpp
#	HeaderFiles/ChoiseDialogPolygon.h
#    SourceFiles/EventFilterManager.cpp
#    SourceFiles/EventFilterManagerGroup.cpp
#    SourceFiles/PolygonTableModel.cpp
#    SourceFiles/PolygonListItem.cpp
#	SourceFiles/MinkowskiSumCalculator.cpp
#    SourceFiles/GuiDrawPolygon.cpp
#    HeaderFiles/EventFilterManager.h
#    HeaderFiles/EventFilterManagerGroup.h
#	HeaderFiles/PolygonTableModel.h#
#	HeaderFiles/PolygonListItem.h
#	HeaderFiles/MinkowskiSumCalculator.h
#	HeaderFiles/MinkowskiSumModified.h
#	HeaderFiles/PolygonItem.h
#	HeaderFiles/QT5/BezierCurves.h
#	HeaderFiles/QT5/GraphicsViewBezierPolygonInput.h
