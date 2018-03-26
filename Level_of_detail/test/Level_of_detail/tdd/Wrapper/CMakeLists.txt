# project(lod)
# cmake_minimum_required(VERSION 2.6)

# find_package(CGAL REQUIRED)
# set(CMAKE_CXX_STANDARD 14)

# if (CGAL_FOUND)
	
# 	include(${CGAL_USE_FILE})
# 	include_directories(BEFORE "../../../../include" "../../../../../../extra-code/libs/OOQP-master/include")

# 	set(CMAKE_CXX_FLAGS "${CMAXE_CXX_FLAGS} -pedantic -Wall -Wextra -Weffc++ -Wno-vla-extension -Wc++11-extensions -Werror -Wno-macro-redefined -Wno-deprecated-declarations")

# 	set(sources main.cpp)

# 	add_library(ma27 STATIC IMPORTED)
# 	set_property(TARGET ma27 PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/../../../../../../extra-code/libs/ma27/src/libma27.a)

# 	add_library(qpbase STATIC IMPORTED)
# 	set_property(TARGET qpbase PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/../../../../../../extra-code/libs/OOQP-master/lib/libooqpbase.a)

# 	add_library(qpdense STATIC IMPORTED)
# 	set_property(TARGET qpdense PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/../../../../../../extra-code/libs/OOQP-master/lib/libooqpdense.a)

# 	add_library(qpgensparse STATIC IMPORTED)
# 	set_property(TARGET qpgensparse PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/../../../../../../extra-code/libs/OOQP-master/lib/libooqpgensparse.a)

# 	add_library(qpgondzio STATIC IMPORTED)
# 	set_property(TARGET qpgondzio PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/../../../../../../extra-code/libs/OOQP-master/lib/libooqpgondzio.a)

# 	add_library(qpsparse STATIC IMPORTED)
# 	set_property(TARGET qpsparse PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/../../../../../../extra-code/libs/OOQP-master/lib/libooqpsparse.a)

# 	add_executable(lod ${sources})
# 	target_link_libraries(lod -L/usr/local/Cellar/gcc/7.2.0/lib/gcc/7/ -lgfortran 
# 	-L/usr/local/Cellar/openblas/0.2.20_1/lib/ -lblas -llapack ma27 qpbase qpdense qpgensparse qpgondzio qpsparse)

# else()
# 	message(STATUS "This program requires the CGAL library and will not be compiled!")
# endif()

# VERSION WITH OPEN CV AND GDAL ---------->

project(lod)
cmake_minimum_required(VERSION 2.6)

# find_package(Eigen3 REQUIRED)
find_package(OpenCV REQUIRED)
find_package(GDAL REQUIRED)
find_package(Boost REQUIRED COMPONENTS filesystem)
find_package(CGAL REQUIRED)

set(CMAKE_CXX_STANDARD 11)

if (CGAL_FOUND)

	include(${CGAL_USE_FILE})
	# include(${EIGEN3_USE_FILE})

	include_directories(BEFORE "../../../../include" ${OpenCV_INCLUDE_DIRS} ${GDAL_INCLUDE_DIR} ${Boost_INCLUDE_DIRS} "../../../../../../extra-code/libs/OOQP-master/include")
	set(CMAKE_CXX_FLAGS "${CMAXE_CXX_FLAGS} -Wno-macro-redefined -Wno-switch -Wno-logical-op-parentheses -Wno-deprecated-declarations")

	set(sources main.cpp
	../../../../include/CGAL/Regularizer/jean_philippe/indexed_event.cpp 
	../../../../include/CGAL/Regularizer/jean_philippe/kinetic_model.cpp 
	../../../../include/CGAL/Regularizer/jean_philippe/line_item.cpp 
	../../../../include/CGAL/Regularizer/jean_philippe/parameters.cpp 
	../../../../include/CGAL/Regularizer/jean_philippe/partition_elements.cpp 
	../../../../include/CGAL/Regularizer/jean_philippe/partition.cpp 
	../../../../include/CGAL/Regularizer/jean_philippe/quadtree_point.cpp 
	../../../../include/CGAL/Regularizer/jean_philippe/segment_ray.cpp 
	../../../../include/CGAL/Regularizer/jean_philippe/segment_tree.cpp 
	../../../../include/CGAL/Regularizer/jean_philippe/geometry.cpp 
	../../../../include/CGAL/Regularizer/jean_philippe/means.cpp 
	../../../../include/CGAL/Regularizer/jean_philippe/regularization_angles.cpp 
	../../../../include/CGAL/Regularizer/jean_philippe/regularization_ordinates.cpp 
	../../../../include/CGAL/Regularizer/jean_philippe/quadratic_problem.cpp 
	../../../../include/CGAL/Regularizer/jean_philippe/regularization_angles_quadratic.cpp 
	../../../../include/CGAL/Regularizer/jean_philippe/regularization_ordinates_quadratic.cpp
	../../../../include/CGAL/Regularizer/jean_philippe/propagation.cpp)

	add_library(ma27 STATIC IMPORTED)
	set_property(TARGET ma27 PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/../../../../../../extra-code/libs/ma27/src/libma27.a)

	add_library(qpbase STATIC IMPORTED)
	set_property(TARGET qpbase PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/../../../../../../extra-code/libs/OOQP-master/lib/libooqpbase.a)

	add_library(qpdense STATIC IMPORTED)
	set_property(TARGET qpdense PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/../../../../../../extra-code/libs/OOQP-master/lib/libooqpdense.a)

	add_library(qpgensparse STATIC IMPORTED)
	set_property(TARGET qpgensparse PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/../../../../../../extra-code/libs/OOQP-master/lib/libooqpgensparse.a)

	add_library(qpgondzio STATIC IMPORTED)
	set_property(TARGET qpgondzio PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/../../../../../../extra-code/libs/OOQP-master/lib/libooqpgondzio.a)

	add_library(qpsparse STATIC IMPORTED)
	set_property(TARGET qpsparse PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/../../../../../../extra-code/libs/OOQP-master/lib/libooqpsparse.a)

	add_executable(lod ${sources})

	target_link_libraries(lod ${OpenCV_LIBS} ${GDAL_LIBRARY} ${Boost_LIBRARIES} -L/usr/local/Cellar/gcc/7.2.0/lib/gcc/7/ -lgfortran 
	-L/usr/local/Cellar/openblas/0.2.20_1/lib/ -lblas -llapack ma27 qpbase qpdense qpgensparse qpgondzio qpsparse)

else()
	message(STATUS "This program requires the CGAL library and will not be compiled!")
endif()