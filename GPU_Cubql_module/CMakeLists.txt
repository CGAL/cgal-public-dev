cmake_minimum_required(VERSION 3.18...3.31) 
project(meshIntersectionRegistration LANGUAGES CXX CUDA)

# ====================================================================
# GLOBAL OPTIONS & BUILD CONFIGURATION
# ====================================================================

option(BUILD_ARCHIVED_TARGETS "Build archived legacy executables" OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL "Export compile commands for IDEs" FORCE)

# Global language standards
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# ====================================================================
# DEPENDENCIES & PACKAGE LOCATIONS
# ====================================================================

find_package(CGAL REQUIRED PATHS "${CMAKE_CURRENT_SOURCE_DIR}/../../cgal-build" NO_DEFAULT_PATH)
find_package(CUDAToolkit REQUIRED)
find_package(Threads REQUIRED)
find_package(TBB REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(READLINE REQUIRED readline)

# Path variables for cuBQL and CUDA context
set(CUBQL_ROOT "/home/yury/Projects/CPPProj/Cubql" CACHE PATH "Path to cuBQL source root")
set(CUBQL_BUILD "${CUBQL_ROOT}/build")

# Shared Link Libraries List
set(SHARED_LINK_LIBS
    CGAL::CGAL
    gmp
    mpfr
    Threads::Threads
    TBB::tbb
    "${CUBQL_BUILD}/libcuBQL_cuda_float2.so"
    "${CUBQL_BUILD}/libcuBQL_cuda_float3.so"
)

# Fetch Polyscope
include(FetchContent)
FetchContent_Declare(
    polyscope
    GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
    GIT_TAG        master
)
FetchContent_MakeAvailable(polyscope)

# ====================================================================
# HELPER FUNCTION
# ====================================================================

function(configure_geometry_target TARGET_NAME)
    set_target_properties(${TARGET_NAME} PROPERTIES
        CUDA_ARCHITECTURES "120"
        CUDA_SEPARABLE_COMPILATION ON
        BUILD_WITH_INSTALL_RPATH TRUE
        INSTALL_RPATH "${CUBQL_BUILD}"
    )

    target_compile_definitions(${TARGET_NAME} PRIVATE CGAL_LINKED_WITH_TBB)

    target_compile_options(${TARGET_NAME} PRIVATE 
        $<$<COMPILE_LANGUAGE:CXX>:-O3 -march=native -DNDEBUG>
        $<$<COMPILE_LANGUAGE:CUDA>:-O3 -DNDEBUG --use_fast_math>
    )

    target_include_directories(${TARGET_NAME} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CUBQL_ROOT}
        ${CUBQL_ROOT}/samples/common
        ${CUDAToolkit_INCLUDE_DIRS}
    )

    target_link_libraries(${TARGET_NAME} PRIVATE 
        CUDA::cudart 
        ${SHARED_LINK_LIBS}
    )
endfunction()

# ====================================================================
# PRIMARY TARGET: MainApp
# ====================================================================

set_source_files_properties(
    apps/mainApplication.cpp 
    src/CPU/TestSuite.cpp 
    PROPERTIES LANGUAGE CUDA
)

add_executable(MainApp 
    apps/mainApplication.cpp 

    src/CPU/TestSuite.cpp
    # Top-Level Entry Points
    src/GPUIntersector/KernelBVHController.cu
    src/GPUIntersector/StandaloneBVHPipeline.cu

    # Traversal & Cross-Checking
    src/traversal/crossCheck.cu
    src/traversal/rapidDescendKernel.cu
    src/traversal/DualTreeStep.cu
    src/traversal/prune_pipeline.cu

    # Single-Tree Traversal & Predicate Dispatcher
    src/single_tree/SingleTreeBatchCommon.cu
    src/single_tree/SingleTreeBatchIntersector.cu

    # Core Infrastructure & Shared Utilities
    src/common/utils.cu
    src/Warmup/cuda_warmup.cu

    # Exact Geometric Predicates
    src/GPU_exact_predicates/TwoSumDoublePredicate.cu
    src/GPU_exact_predicates/InexactBigIntFilter.cu
    src/GPU_exact_predicates/ExactBigIntPredicate.cu
    src/GPU_exact_predicates/OriginShiftedFloatPredicate.cu
    src/GPU_exact_predicates/DoublePredicate.cu
)

configure_geometry_target(MainApp)

# Target-specific additions for MainApp
target_include_directories(MainApp PRIVATE ${READLINE_INCLUDE_DIRS})
target_link_libraries(MainApp PRIVATE ${READLINE_LIBRARIES} polyscope)

# ====================================================================
# ARCHIVED / LEGACY TARGETS (OPTIONAL BUILD)
# ====================================================================

if(BUILD_ARCHIVED_TARGETS)
    message(STATUS "Building archived/legacy targets enabled.")

    # TARGET: meshIntersectionRegistration
    add_executable(meshIntersectionRegistration apps/archive/main.cpp src/registration/kernels.cu)
    configure_geometry_target(meshIntersectionRegistration)

    # TARGET: meshIntersectionBatch
    add_executable(meshIntersectionBatch apps/archive/main_massive.cpp src/registration/kernels.cu)
    configure_geometry_target(meshIntersectionBatch)

    # TARGET: meshIntersectionCustom
    add_executable(meshIntersectionCustom apps/archive/mainV2.cpp src/GPU_exact_predicates/kernelsV2.cu src/GPU_exact_predicates/InexactBigIntFilter.cu)
    configure_geometry_target(meshIntersectionCustom)

    # TARGET: meshVoxelizer
    add_executable(meshVoxelizer apps/archive/mainVoxel.cpp src/voxelizer/kernelsVoxel.cu)
    configure_geometry_target(meshVoxelizer)

    # TARGET: testBVH
    add_executable(testBVH apps/archive/mainTestBVH.cpp src/archive/kernelsTestBVH.cu src/traversal/rapidDescendKernel.cu src/traversal/crossCheck.cu src/archive/volumeSanityCheck.cu src/archive/batchedCrossIntersection.cu src/archive/batchedCrossIntersectionBruteV2.cu src/traversal/DualTreeStep.cu src/traversal/prune_pipeline.cu src/GPU_exact_predicates/InexactBigIntFilter.cu)
    configure_geometry_target(testBVH)

    # TARGET: dualTreeExecution
    add_executable(dualTreeExecution apps/archive/mainDualTreeTestLauncher.cpp src/traversal/crossCheck.cu src/GPU_exact_predicates/InexactBigIntFilter.cu src/DualTreeTraversal/DualTreeLauncher.cu src/DualTreeTraversal/DualTreeMeshIntersection.cu)
    configure_geometry_target(dualTreeExecution)

    # TARGET: testBVHV2
    add_executable(testBVHV2 apps/archive/mainTestBVHV2.cpp src/archive/kernelsTestBVHAlternative.cu src/traversal/rapidDescendKernel.cu src/archive/batchedCrossIntersectionV2.cu src/archive/batchedCrossIntersectionV3.cu src/single_tree/SingleTreeBatchCommon.cu src/traversal/DualTreeStep.cu src/traversal/prune_pipeline.cu src/Warmup/cuda_warmup.cu src/GPU_exact_predicates/ExactBigIntPredicate.cu src/common/utils.cu)
    configure_geometry_target(testBVHV2)

endif()