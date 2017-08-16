project(mylogtest)

cmake_minimum_required(VERSION 2.6)

find_package(CGAL REQUIRED)

set(CMAKE_CXX_STANDARD 11)

if (CGAL_FOUND) 

	include(${CGAL_USE_FILE})

	include_directories(BEFORE "../../../../include")

	include_directories($ENV{GMOCK_HOME}/googlemock/include $ENV{GMOCK_HOME}/googletest/include)

	link_directories($ENV{GMOCK_HOME}/googlemock/build $ENV{GMOCK_HOME}/googletest/build)

	set(CMAKE_CXX_FLAGS "${CMAXE_CXX_FLAGS} -Wall -Wc++11-extensions")

	set(sources main.cpp MyLogTest.cpp)

	add_executable(mylogtest ${sources})

	target_link_libraries(mylogtest pthread)
	target_link_libraries(mylogtest gmock)
	target_link_libraries(mylogtest gtest)

else()

	message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()