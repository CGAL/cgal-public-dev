
##############################################################################
# There are two ways to add include directories to the NVCC command
# line:

cmake_minimum_required(VERSION 2.8)

project(gpu_symbolic)

# whether to compile for 32- or 64-bit architecture (default: 0)
set(CUMP_USE_32_BITS 0)
# whether to use 32-bit or 24-bit moduli set (default: 1)
set(USE_32_BIT_MODULI_SET 0)

if(CUMP_USE_32_BITS)

endif()

set(CUDA_TOOLKIT_ROOT_DIR "$ENV{CUDA_TOOLKIT_ROOT_DIR}")
set(COMMON_INC_DIR "$ENV{HOME}/work/GPU_programs/common")
set(GCD_DIR "$ENV{HOME}/work/GPU_programs/gpu_gcd")
set(RES_DIR "$ENV{HOME}/work/GPU_programs/gpu_resultant")

if(USE_32_BIT_MODULI_SET)

set(CUDA_ARCH_FLAG "sm_20")
set(MODULI_SET_FLAG "-DCUMP_USE_32BIT_MODULI_SET=1")

else()

set(CUDA_ARCH_FLAG "sm_13")
set(MODULI_SET_FLAG "-DCUMP_USE_32BIT_MODULI_SET=0")

endif()

find_package(CUDA QUIET REQUIRED)
find_package(CGAL QUIET COMPONENTS Core)

include("cgal_cmake_lib_defs")

set(CMAKE_BUILD_TYPE Release)

# NOTE NOTE NOTE: do not forget about '-pg' option !!!!

# -keep --ptxas-options;--maxrregcount=128
# set(CUDA_NVCC_FLAGS "--compiler-options;-mtune=core2;-mfpmath=sse")
set(CUDA_NVCC_FLAGS "--compiler-options;-fno-strict-aliasing;--compiler-options;-fPIC;-arch;${CUDA_ARCH_FLAG};--ptxas-options;-v")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MODULI_SET_FLAG}")

# Set CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE when you want to add the same .cu
# file to multiple targets.
#set(CUDA_PROPAGATE_HOST_FLAGS OFF)
set(CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE OFF)
set(CUDA_VERBOSE_BUILD ON)

if(CUMP_USE_32_BITS)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m32")
set(CUDA_64_BIT_DEVICE_CODE OFF)
set(CUDA_ARCH_TYPE "32") 
message(STATUS "32-bit device code")
else()
set(CUDA_64_BIT_DEVICE_CODE ON)
set(CUDA_ARCH_TYPE "64")
message(STATUS "64-bit device code")
endif()

include_directories( BEFORE ${COMMON_INC_DIR} ${COMMON_INC_DIR}/include ${GCD_DIR} ${GCD_DIR}/include ${RES_DIR} ${RES_DIR}/include include .)

cuda_compile(CUDA_FILES ${GCD_DIR}/gcd_algorithm_gpu.cu
    ${RES_DIR}/resultant_algorithm_gpu.cu device_manager.cu)

add_library(gpu_symbolic_${CUDA_ARCH_FLAG}_${CUDA_ARCH_TYPE} SHARED
    ${CUDA_FILES}
    ${GCD_DIR}/gcd_algorithm_host.C
    ${RES_DIR}/resultant_algorithm_host.C
    ${COMMON_INC_DIR}/modular_arithm.C)

# Add a special target to clean nvcc generated files.
CUDA_BUILD_CLEAN_TARGET()

