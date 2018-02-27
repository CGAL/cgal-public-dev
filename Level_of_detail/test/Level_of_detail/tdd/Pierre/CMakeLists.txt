# project(lod)
# cmake_minimum_required(VERSION 2.6)

# find_package(CGAL REQUIRED)
# set(CMAKE_CXX_STANDARD 11)

# if (CGAL_FOUND)
# 	include(${CGAL_USE_FILE})

# 	include_directories(BEFORE "../../../../include" "C:/dev/") 
	
# 	if (WIN32)
# 		set(CMAKE_CXX_FLAGS "${CMAXE_CXX_FLAGS} /EHsc")
# 	else()
# 		set(CMAKE_CXX_FLAGS "${CMAXE_CXX_FLAGS}")
# 	endif()

# 	set(sources main.cpp)
# 	add_executable(lod ${sources})

# else()
# 	message(STATUS "This program requires the CGAL library and will not be compiled!")
# endif()
