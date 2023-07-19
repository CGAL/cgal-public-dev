@echo off
set folder="build\"
if not exist %folder% (
   echo "no build folder. terminating..."
   exit /b
) 

cmake --build build 2> error.cmake.txt
