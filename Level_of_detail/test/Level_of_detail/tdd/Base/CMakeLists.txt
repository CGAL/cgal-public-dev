project(basetest)

cmake_minimum_required(VERSION 2.6)

find_package(CGAL REQUIRED)

set(CMAKE_CXX_STANDARD 14)

if (CGAL_FOUND) 

	include(${CGAL_USE_FILE})

	include_directories(BEFORE "../../../../include")

	include_directories($ENV{GMOCK_HOME}/googlemock/include $ENV{GMOCK_HOME}/googletest/include)

	link_directories($ENV{GMOCK_HOME}/googlemock/build $ENV{GMOCK_HOME}/googletest/build)

	set(CMAKE_CXX_FLAGS "${CMAXE_CXX_FLAGS} -pedantic -Wall -Wextra -Weffc++ -Wno-vla-extension -Wc++11-extensions -Werror")

	set(sources main.cpp Level_of_detail_base_test.cpp)

	add_executable(basetest ${sources})

	target_link_libraries(basetest pthread)
	target_link_libraries(basetest gmock)
	target_link_libraries(basetest gtest)

else()

	message(STATUS "This program requires the CGAL library and will not be compiled!")

endif()