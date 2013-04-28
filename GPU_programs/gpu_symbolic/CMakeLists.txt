
##############################################################################
# There are two ways to add include directories to the NVCC command
# line:

cmake_minimum_required(VERSION 2.8)

project(gpu_symbolic)

set(CUMP_USE_32_BITS 0)
set(USE_32_BIT_MODULI_SET 1)

set(CUDA_TOOLKIT_ROOT_DIR "$ENV{CUDA_TOOLKIT_ROOT_DIR}")

if(CUMP_USE_32_BITS)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m32 -fPIC")
else()
endif()

set(COMMON_INC_DIR "$ENV{HOME}/work/GPU_programs/common")
set(GPU_SYMB_DIR "$ENV{HOME}/work/GPU_programs/gpu_symbolic")
set(GCD_DIR "$ENV{HOME}/work/GPU_programs/gpu_gcd")
set(RES_DIR "$ENV{HOME}/work/GPU_programs/gpu_resultant")

if(USE_32_BIT_MODULI_SET)

set(CUDA_ARCH_FLAG "sm_20")
set(MODULI_SET_FLAG "-DCUMP_USE_32BIT_MODULI_SET=1")

else()

set(CUDA_ARCH_FLAG "sm_13")
set(MODULI_SET_FLAG "-DCUMP_USE_32BIT_MODULI_SET=0")

endif()

if(CUMP_USE_32_BITS)
set(SONAME "${CUDA_ARCH_FLAG}_32")
else()
set(SONAME "${CUDA_ARCH_FLAG}_64")
endif()

find_package(CUDA QUIET REQUIRED)
find_package(CGAL QUIET COMPONENTS Core)

FIND_PACKAGE( Boost REQUIRED COMPONENTS thread )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

include("cgal_cmake_lib_defs")

set(CMAKE_BUILD_TYPE Release)


# SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE) 

# set(CUDA_NVCC_FLAGS "--compiler-options;-mtune=core2;-mfpmath=sse")
set(CUDA_NVCC_FLAGS "--compiler-options;-fno-strict-aliasing;-arch;${CUDA_ARCH_FLAG};--ptxas-options;-v;--compiler-options;-fPIC")
# set(CUDA_INCLUDE_DIRS  "${CUDA_CUT_INCLUDE_DIR}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MODULI_SET_FLAG}")

# Set CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE when you want to add the same .cu
# file to multiple targets.
#set(CUDA_PROPAGATE_HOST_FLAGS OFF)
set(CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE OFF)
set(CUDA_VERBOSE_BUILD ON)

if(CUMP_USE_32_BITS)
set(CUDA_64_BIT_DEVICE_CODE OFF)
else()
set(CUDA_64_BIT_DEVICE_CODE ON)
endif()


include_directories( BEFORE ${COMMON_INC_DIR} ${COMMON_INC_DIR}/include ${GCD_DIR} ${GCD_DIR}/include ${RES_DIR} ${RES_DIR}/include include .)

cuda_compile(CUDA_FILES ${GCD_DIR}/gcd_algorithm_gpu.cu
    ${RES_DIR}/resultant_algorithm_gpu.cu device_manager.cu)

cuda_add_library(gpu_symbolic_${SONAME} SHARED 
    ${CUDA_FILES}
    ${GCD_DIR}/gcd_algorithm_host.C
    ${RES_DIR}/resultant_algorithm_host.C
    ${COMMON_INC_DIR}/modular_arithm.C)

TARGET_LINK_LIBRARIES(gpu_symbolic_${SONAME}  ${Boost_LIBRARIES} )
   
#    target_link_libraries(gpu_gcd ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${CUDA_CUDART_LIBRARY} -L${CUDA_DRV_DIR} -lcuda -lstdc++ -L/opt/gcc-4.4.7/lib64)

# Add a special target to clean nvcc generated files.
CUDA_BUILD_CLEAN_TARGET()

