cmake_minimum_required(VERSION 3.18...3.31) 
project(meshIntersectionRegistration LANGUAGES CXX CUDA)

# ====================================================================
# DEPENDENCIES & PACKAGE LOCATIONS
# ====================================================================

# 1. Locate the CGAL build configuration from your parent build directory
find_package(CGAL REQUIRED PATHS "${CMAKE_CURRENT_SOURCE_DIR}/../../cgal-build" NO_DEFAULT_PATH)

# 2. Locate the modern CUDA Toolkit wrapper to populate CUDA:: targets
find_package(CUDAToolkit REQUIRED)

# 3. Define path variables for cuBQL and CUDA context
set(CUBQL_ROOT "/home/yury/Projects/CPPProj/Cubql" CACHE PATH "Path to cuBQL source root")
set(CUBQL_BUILD "${CUBQL_ROOT}/build")

# 4. Find threading dependencies
find_package(Threads REQUIRED)
find_package(TBB REQUIRED)

# Shared Link Libraries List to avoid duplication across targets
set(SHARED_LINK_LIBS
    CGAL::CGAL
    gmp
    mpfr
    Threads::Threads
    TBB::tbb
    "${CUBQL_BUILD}/libcuBQL_cuda_float2.so"
    "${CUBQL_BUILD}/libcuBQL_cuda_float3.so"
)

# ====================================================================
# TARGET 1: Your Original Project (PRISTINE: main.cpp + kernels.cu)
# ====================================================================
add_executable(meshIntersectionRegistration apps/main.cpp src/registration/kernels.cu)

set_target_properties(meshIntersectionRegistration PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
    CUDA_STANDARD 17
    CUDA_STANDARD_REQUIRED ON
    CUDA_ARCHITECTURES "120"
    CUDA_SEPARABLE_COMPILATION ON
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH "${CUBQL_BUILD}"
)

target_compile_definitions(meshIntersectionRegistration PRIVATE CGAL_LINKED_WITH_TBB)
target_compile_options(meshIntersectionRegistration PRIVATE 
    $<$<COMPILE_LANGUAGE:CXX>:-O3 -march=native -DNDEBUG>
    $<$<COMPILE_LANGUAGE:CUDA>:-O3 -DNDEBUG --use_fast_math>
)
target_include_directories(meshIntersectionRegistration PRIVATE
    include
    ${CUBQL_ROOT}
    ${CUBQL_ROOT}/samples/common
    ${CUDAToolkit_INCLUDE_DIRS}
)
target_link_libraries(meshIntersectionRegistration PRIVATE ${SHARED_LINK_LIBS})


# ====================================================================
# TARGET 2: Massive Test Runner Variant 
# ====================================================================
add_executable(meshIntersectionBatch apps/main_massive.cpp src/registration/kernels.cu)

set_target_properties(meshIntersectionBatch PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
    CUDA_STANDARD 17
    CUDA_STANDARD_REQUIRED ON
    CUDA_ARCHITECTURES "120"
    CUDA_SEPARABLE_COMPILATION ON
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH "${CUBQL_BUILD}"
)

target_compile_definitions(meshIntersectionBatch PRIVATE CGAL_LINKED_WITH_TBB)
target_compile_options(meshIntersectionBatch PRIVATE 
    $<$<COMPILE_LANGUAGE:CXX>:-O3 -march=native -DNDEBUG>
    $<$<COMPILE_LANGUAGE:CUDA>:-O3 -DNDEBUG --use_fast_math>
)
target_include_directories(meshIntersectionBatch PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    include
    ${CUBQL_ROOT}
    ${CUBQL_ROOT}/samples/common
    ${CUDAToolkit_INCLUDE_DIRS}
)
target_link_libraries(meshIntersectionBatch PRIVATE ${SHARED_LINK_LIBS})


# ====================================================================
# TARGET 3: Quick Precision Pair Checker (pair_check.cpp)
# ====================================================================
add_executable(pairCheck apps/pair_check.cpp)

set_target_properties(pairCheck PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
)
target_link_libraries(pairCheck PRIVATE CGAL::CGAL gmp mpfr)


# ====================================================================
# TARGET 4: New Partitioned Dual-List Pipeline Sandbox (YOUR MAIN PROJECT)
# ====================================================================
add_executable(meshIntersectionCustom apps/mainV2.cpp src/custom_pipeline/kernelsV2.cu)

set_target_properties(meshIntersectionCustom PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
    CUDA_STANDARD 17
    CUDA_STANDARD_REQUIRED ON
    CUDA_ARCHITECTURES "120"
    CUDA_SEPARABLE_COMPILATION ON
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH "${CUBQL_BUILD}"
)

target_compile_definitions(meshIntersectionCustom PRIVATE CGAL_LINKED_WITH_TBB)

# REMOVED --use_fast_math HERE TO PRESERVE STRICT IEEE-754 DIRECTED ROUNDING FOR INTERVALS
target_compile_options(meshIntersectionCustom PRIVATE 
    $<$<COMPILE_LANGUAGE:CXX>:-O3 -march=native -DNDEBUG>
    $<$<COMPILE_LANGUAGE:CUDA>:-O3 -DNDEBUG>
)
target_include_directories(meshIntersectionCustom PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    include
    ${CUBQL_ROOT}
    ${CUBQL_ROOT}/samples/common
)
target_link_libraries(meshIntersectionCustom PRIVATE CUDA::cudart ${SHARED_LINK_LIBS})


# ====================================================================
# TARGET 5: Dual-Pass Voxel Space Partitioning Engine
# ====================================================================
add_executable(meshVoxelizer apps/mainVoxel.cpp src/voxelizer/kernelsVoxel.cu)

set_target_properties(meshVoxelizer PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
    CUDA_STANDARD 17
    CUDA_STANDARD_REQUIRED ON
    CUDA_ARCHITECTURES "120"
    CUDA_SEPARABLE_COMPILATION ON
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH "${CUBQL_BUILD}"
)

target_compile_definitions(meshVoxelizer PRIVATE CGAL_LINKED_WITH_TBB)
target_compile_options(meshVoxelizer PRIVATE 
    $<$<COMPILE_LANGUAGE:CXX>:-O3 -march=native -DNDEBUG>
    $<$<COMPILE_LANGUAGE:CUDA>:-O3 -DNDEBUG --use_fast_math>
)
target_include_directories(meshVoxelizer PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    include
    ${CUBQL_ROOT}
    ${CUBQL_ROOT}/samples/common
)
target_link_libraries(meshVoxelizer PRIVATE CUDA::cudart ${SHARED_LINK_LIBS})


# ====================================================================
# TARGET 6: TBB CPU-Only Spatial Voxelization SandBox
# ====================================================================
add_executable(meshVoxelizerCPU src/voxelizer/mainCPUVoxel.cpp)

set_target_properties(meshVoxelizerCPU PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
)
target_compile_options(meshVoxelizerCPU PRIVATE -O3 -march=native -DNDEBUG)
target_include_directories(meshVoxelizerCPU PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} include)
target_link_libraries(meshVoxelizerCPU PRIVATE Threads::Threads TBB::tbb)


# ====================================================================
# TARGET 7: Single-Mesh BVH Tree Layout Diagnostic Tester
# ====================================================================
add_executable(testBVH 
    apps/mainTestBVH.cpp 
    src/testBVH/kernelsTestBVH.cu 
    src/testBVH/rapidDescendKernel.cu
    src/testBVH/crossCheck.cu
    src/testBVH/volumeSanityCheck.cu
    src/testBVH/batchedCrossIntersection.cu       
)

set_target_properties(testBVH PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
    CUDA_STANDARD 17
    CUDA_STANDARD_REQUIRED ON
    CUDA_ARCHITECTURES "120"
    CUDA_SEPARABLE_COMPILATION ON
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH "${CUBQL_BUILD}"
)

target_compile_definitions(testBVH PRIVATE CGAL_LINKED_WITH_TBB)

# Kept --use_fast_math enabled here since we aren't using floating-point interval arithmetic 
target_compile_options(testBVH PRIVATE 
    $<$<COMPILE_LANGUAGE:CXX>:-O3 -march=native -DNDEBUG>
    $<$<COMPILE_LANGUAGE:CUDA>:-O3 -DNDEBUG --use_fast_math>
)

target_include_directories(testBVH PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/GPU_Cubql_module          # Path tracking to your root module
    ${CMAKE_CURRENT_SOURCE_DIR}/GPU_Cubql_module/include  # Path tracking to where sm_builder_v2.h sits
    include
    ${CUBQL_ROOT}
    ${CUBQL_ROOT}/samples/common
    ${CUDAToolkit_INCLUDE_DIRS}
)

target_link_libraries(testBVH PRIVATE CUDA::cudart ${SHARED_LINK_LIBS})