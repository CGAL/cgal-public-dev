
##############################################################################
# There are two ways to add include directories to the NVCC command
# line:

cmake_minimum_required(VERSION 2.8)

project(test_CUDA)

# INCLUDE (CheckCXXSourceCompiles)
INCLUDE (CheckCXXSourceRuns)
# INCLUDE (CheckCSourceCompiles)

# set(CUDA_TOOLKIT_ROOT_DIR "$ENV{HOME}/cuda1.0/cuda")

find_package(CUDA QUIET REQUIRED)
set(CMAKE_BUILD_TYPE Release)

set(CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE OFF)
set(CUDA_VERBOSE_BUILD ON)
set(CUDA_64_BIT_DEVICE_CODE OFF)

# set(src /home/emeliyan/work/GPU_programs/test_CUDA/test_CUDA.C)
# message(STATUS "cuda inc =  ${CUDA_INCLUDE_DIRS}")

# CMAKE_REQUIRED_FLAGS = string of compile command line flags
#   CMAKE_REQUIRED_DEFINITIONS = list of macros to define (-DFOO=bar)
#   CMAKE_REQUIRED_INCLUDES = list of include directories
#   CMAKE_REQUIRED_LIBRARIES = list of libraries to link

set(CMAKE_REQUIRED_INCLUDES ${CUDA_INCLUDE_DIRS})
set(CMAKE_REQUIRED_LIBRARIES ${CUDA_LIBRARIES})
FILE(READ test_CUDA.C src)

# set(f "this; ...is: a; test")
# string(REPLACE E F f "${src}")

#   string(REPLACE <match_string>
#          <replace_string> <output variable>
#          <input> [<input>...])
# 
# CHECK_CXX_SOURCE_COMPILES("${src}" DDD)

CHECK_CXX_SOURCE_RUNS("${src}" CUDA_CHECK)

if(${CUDA_CHECK_COMPILED})

    message(STATUS  "test_CUDA compiled successfully")
    MATH(EXPR CUDA_MAJOR '${CUDA_CHECK_EXITCODE}>>4')
    MATH(EXPR CUDA_MINOR '${CUDA_CHECK_EXITCODE}%16')

    if(${CUDA_MAJOR} EQUAL 0)
        message(STATUS "\nThere is not device supporting CUDA\n")        
    else()
        message(STATUS  "\n========= CUDA version found: ${CUDA_MAJOR}.${CUDA_MINOR} ============\n")
    endif()

else()
    message(STATUS  "test_CUDA failed: ${OUTPUT}")
endif()

 # CUDA_ADD_EXECUTABLE(test_CUDA
#     test_CUDA.C)

# Add a special target to clean nvcc generated files.
CUDA_BUILD_CLEAN_TARGET()

