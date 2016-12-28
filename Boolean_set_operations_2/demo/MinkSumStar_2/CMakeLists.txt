cmake_minimum_required(VERSION 2.8)

project(MinkowskiSumStar_demo)

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

qt5_add_resources(UI_RESOURCES
    ResourceFiles/icons.qrc)

qt5_wrap_ui(UI_HEADERS
    UiFiles/MainWindow.ui
	UiFiles/ChoiseDialogPolygon.ui)

add_executable(MinkowskiSumStar_demo
	Boolean_set_operations_2.cpp
    SourceFiles/MainWindow.cpp
    SourceFiles/EventFilterManager.cpp
    SourceFiles/EventFilterManagerGroup.cpp
    SourceFiles/PolygonTableModel.cpp
    SourceFiles/PolygonListItem.cpp
	SourceFiles/MinkowskiSumCalculator.cpp
	SourceFiles/ChoiseDialogPolygon.cpp
    HeaderFiles/MainWindow.h
    HeaderFiles/EventFilterManager.h
    HeaderFiles/EventFilterManagerGroup.h
	HeaderFiles/PolygonTableModel.h
	HeaderFiles/PolygonListItem.h
	HeaderFiles/MinkowskiSumCalculator.h
	HeaderFiles/MinkowskiSumModified.h
	HeaderFiles/ChoiseDialogPolygon.h
	HeaderFiles/PolygonItem.h
	HeaderFiles/Typedefs.h
    ${UI_HEADERS}
    ${UI_RESOURCES})

# Link with Qt libraries
target_link_libraries(MinkowskiSumStar_demo ${QT_LIBRARIES})

# Link with CGAL
target_link_libraries(MinkowskiSumStar_demo ${CGAL_LIBRARIES}  ${CGAL_3RD_PARTY_LIBRARIES})

qt5_use_modules(MinkowskiSumStar_demo
    Widgets
    Script
    Svg
	OpenGL)

#set_target_properties(MinkowskiSumStar_demo PROPERTIES LINK_FLAGS "/SUBSYSTEM:WINDOWS /ENTRY:mainCRTStartup")
