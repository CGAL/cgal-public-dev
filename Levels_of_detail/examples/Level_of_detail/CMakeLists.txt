project(Level_of_detail_example)

cmake_minimum_required(VERSION 2.6)
find_package(CGAL REQUIRED)

set(CMAKE_CXX_STANDARD 11)
if (CGAL_FOUND)

	include(${CGAL_USE_FILE})
	include_directories(BEFORE "../../include")

	set(CMAKE_CXX_FLAGS "${CMAXE_CXX_FLAGS} -Wno-deprecated-declarations -Wno-user-defined-warnings -Werror")
	set(sources Level_of_detail_example.cpp)

	add_executable(Level_of_detail_example ${sources})
	target_link_libraries(Level_of_detail_example)

else()
	message(STATUS "This program requires the CGAL library and will not be compiled!")
endif()