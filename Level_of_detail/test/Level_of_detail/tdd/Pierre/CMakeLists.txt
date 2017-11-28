project(lod)
cmake_minimum_required(VERSION 2.6)

find_package(CGAL REQUIRED)
set(CMAKE_CXX_STANDARD 11)

if (CGAL_FOUND)
	include(${CGAL_USE_FILE})

	include_directories(BEFORE "../../../../include" "C:/dev/") 
	
	# windows
	# set(CMAKE_CXX_FLAGS "${CMAXE_CXX_FLAGS} /EHsc")
	
	# mac os
	set(CMAKE_CXX_FLAGS "${CMAXE_CXX_FLAGS}")

	set(sources main.cpp)
	add_executable(lod ${sources})

else()
	message(STATUS "This program requires the CGAL library and will not be compiled!")
endif()