@echo off

if not defined CGAL_DIR (
   echo "You must define the CGAL_DIR environment variable. Terminating..."
   exit /b
)

if not defined Boost_INCLUDE_DIR (
   echo "You must define the Boost_INCLUDE_DIR environment variable. Terminating..."   
   exit /b
)

if not defined CMAKE_TOOLCHAIN_FILE (
   echo "You must define the CMAKE_TOOLCHAIN_FILE environment variable. Terminating..."
   exit /b
)

set folder="build\"
if not exist %folder% (
   echo "no build folder. creating..."
   mkdir build
)

set folder="src\"
if not exist %folder% (
   echo "no src folder. terminating..."
   exit /b
)

cmake -B .\build -S .\src --preset=default