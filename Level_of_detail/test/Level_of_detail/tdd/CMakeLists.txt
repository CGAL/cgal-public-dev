project(logtest)

cmake_minimum_required(VERSION 2.6)

find_package(CGAL REQUIRED)

if (CGAL_FOUND) 

	include( ${CGAL_USE_FILE} )

	include_directories (BEFORE "../../../include")

	include_directories($ENV{GMOCK_HOME}/googlemock/include $ENV{GMOCK_HOME}/googletest/include)

	link_directories($ENV{GMOCK_HOME}/googlemock/build $ENV{GMOCK_HOME}/googletest/build)

	set(CMAKE_CXX_FLAGS "${CMAXE_CXX_FLAGS} -Wall")

	set(sources main.cpp LogTest.cpp)

	add_executable(logtest ${sources})

	target_link_libraries(logtest pthread)
	target_link_libraries(logtest gmock)
	target_link_libraries(logtest gtest)

else()

	message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()