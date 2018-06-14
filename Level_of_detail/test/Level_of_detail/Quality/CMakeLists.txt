project(quality)
cmake_minimum_required(VERSION 2.6)

find_package(Eigen3 REQUIRED)
find_package(OpenCV REQUIRED)
find_package(GDAL REQUIRED)
find_package(Boost REQUIRED COMPONENTS filesystem)
find_package(CGAL REQUIRED)

set(CMAKE_CXX_STANDARD 11)
if (CGAL_FOUND)

	include(${CGAL_USE_FILE})
	include_directories(BEFORE 
	"../../../include"
	${EIGEN3_INCLUDE_DIR} 
	${OpenCV_INCLUDE_DIRS} 
	${GDAL_INCLUDE_DIR} 
	${Boost_INCLUDE_DIRS} 
	"../../../../../libs/ooqp/include" 
	"../../../../../libs/kinetic/include")

	set(CMAKE_CXX_FLAGS "${CMAXE_CXX_FLAGS} -Wno-deprecated-declarations -Wno-user-defined-warnings -Werror")
	set(sources main.cpp)

	add_library(libkinetic STATIC IMPORTED)
	set_property(TARGET libkinetic PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/../../../../../libs/kinetic/build-release/libkinetic.a)

	add_library(ma27 STATIC IMPORTED)
	set_property(TARGET ma27 PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/../../../../../libs/ma27/src/libma27.a)

	add_library(qpbase STATIC IMPORTED)
	set_property(TARGET qpbase PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/../../../../../libs/ooqp/lib/libooqpbase.a)

	add_library(qpdense STATIC IMPORTED)
	set_property(TARGET qpdense PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/../../../../../libs/ooqp/lib/libooqpdense.a)

	add_library(qpgensparse STATIC IMPORTED)
	set_property(TARGET qpgensparse PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/../../../../../libs/ooqp/lib/libooqpgensparse.a)

	add_library(qpgondzio STATIC IMPORTED)
	set_property(TARGET qpgondzio PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/../../../../../libs/ooqp/lib/libooqpgondzio.a)

	add_library(qpsparse STATIC IMPORTED)
	set_property(TARGET qpsparse PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/../../../../../libs/ooqp/lib/libooqpsparse.a)

	add_executable(quality ${sources})

	target_link_libraries(quality ${OpenCV_LIBS} ${GDAL_LIBRARY} ${Boost_LIBRARIES} -L/usr/local/Cellar/gcc/8.1.0/lib/gcc/8/ -lgfortran 
	-L/usr/local/Cellar/openblas/0.2.20_1/lib/ -lblas -llapack ma27 qpbase qpdense qpgensparse qpgondzio qpsparse libkinetic)

else()
	message(STATUS "This program requires the CGAL library and will not be compiled!")
endif()