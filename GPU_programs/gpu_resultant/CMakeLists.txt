
cmake_minimum_required(VERSION 2.8)

project(gpu_resultant)

set(CUMP_USE_32_BITS 0)

set(CUDA_TOOLKIT_ROOT_DIR "$ENV{CUDA_TOOLKIT_ROOT_DIR}")
set(CUDA_DRV_DIR "/usr/lib64/nvidia")

if(CUMP_USE_32_BITS)
set(CMAKE_CXX_FLAGS "${CMAKE_XXC_FLAGS} -m32")
else()
endif()

set(COMMON_INC_DIR "$ENV{HOME}/work/GPU_programs/common")
set(GPU_SYMB_DIR "$ENV{HOME}/work/GPU_programs/gpu_symbolic")

find_package(CUDA QUIET REQUIRED)
find_package(CGAL QUIET COMPONENTS Core)

if(NOT CUMP_COMPILE_LIBRARY) 

include( ${CGAL_USE_FILE} )

else()

include("cgal_cmake_lib_defs")

endif()

# set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_DEBUG "-frounding-math")

set(CUDA_ARCH_FLAG "sm_20")

# set(CUDA_NVCC_FLAGS "--compiler-options;-mtune=core2;-mfpmath=sse")
set(CUDA_NVCC_FLAGS "--compiler-options;-fno-strict-aliasing;-arch;${CUDA_ARCH_FLAG};--ptxas-options;-v")
# set(CUDA_INCLUDE_DIRS  "${CUDA_CUT_INCLUDE_DIR}")

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

include_directories( BEFORE $ENV{HOME}/work/CGAL/branches/experimental-packages/Algebraic_kernel_d/include ${GPU_SYMB_DIR} ${COMMON_INC_DIR} ${COMMON_INC_DIR}/include include .)

cuda_compile(CUDA_FILES resultant_algorithm_gpu.cu
         ${GPU_SYMB_DIR}/device_manager.cu)

if(CUMP_COMPILE_LIBRARY) 

add_library(gpu_resultant_${CUDA_ARCH_FLAG}
    ${CUDA_FILES}
    resultant_algorithm_host.C
    ${COMMON_INC_DIR}/modular_arithm.C)
else()

add_executable(gpu_resultant
    main.C
    ${CUDA_FILES}
    resultant_algorithm_host.C
    ${COMMON_INC_DIR}/modular_arithm.C)

# set(CUDA_CUDART_LIBRARY "libcudart.so")

target_link_libraries(gpu_resultant ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${CUDA_CUDART_LIBRARY} -lcuda -lstdc++)


endif()

# target_link_libraries(gpu_resultant ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES}
# ${CUDA_CUT_LIBRARY} ${CUDA_CUDART_LIBRARY}
# )

# Add a special target to clean nvcc generated files.
CUDA_BUILD_CLEAN_TARGET()

